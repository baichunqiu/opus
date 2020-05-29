#include "../include/beamformer.h"
#include "../include/complex.h"
#include "../include/icm_fft.h"
#include "../include/ring_buffer.h"
#include "../include/array_utils.h"

#define kFftSize                    512
#define kFreqBins                   (kFftSize/2+1)
#define kLowMeanStartHz             200
#define kLowMeanEndHz               400
#define kSpeedOfSoundMeterSeconds   343.0
// 双麦克参数
#define kBalance                    0.4
#define kBeamwidthConstant          0.01
// 4麦参数
//#define kBalance                    0.9
//#define kBeamwidthConstant          0.000002
#define kMaskMinimum                0.01
#define kMaskTimeSmoothAlpha        0.2
#define kMaskFrequencySmoothAlpha   0.6
#define S16_MAX                     32768.0

#define TAG "log" // 这个是自定义的LOG的标识
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG,TAG ,__VA_ARGS__) // 定义LOGD类型
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG ,__VA_ARGS__) // 定义LOGI类型
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN,TAG ,__VA_ARGS__) // 定义LOGW类型
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR,TAG ,__VA_ARGS__) // 定义LOGE类型
#define LOGF(...) __android_log_print(ANDROID_LOG_FATAL,TAG ,__VA_ARGS__) // 定义LOGF类型

struct Beamformer
{
    OsInt32 fs;
    OsInt32 frmsize;
    OsFloat targetAngle;
    OsFloat interfAngle;
    OsInt32 low_mean_start_bin;
    OsInt32 low_mean_end_bin;
    OsInt32 high_mean_start_bin;
    OsInt32 high_mean_end_bin;
    OsInt32 num_output_channels;
    SensorPoint geometry[kNumMics];
    OsInt32 init_delay;
    OsInt32 frame_offset;
    OsInt32 shift_amount;

    OsFloat     *window;
    OsFloat     *x_buf[kNumMics];
    OsFloat     *y_buf;
    OsFloat     *norm_in_buf[kNumMics];
    OsFloat     *norm_out_buf;
    RingBuffer  *in_buf[kNumMics];
    OsFloat     *out_buf;
    complex     *in_block[kNumMics];
    complex     *out_block;
    fft_t       *fft;

    complex delay_sum_masks[kFreqBins][kNumMics];
    complex normalized_delay_sum_masks[kFreqBins][kNumMics];
    complex target_cov_mats[kFreqBins][kNumMics][kNumMics];
    complex interf_cov_mats[kFreqBins][kNumMics][kNumMics];
    complex reflected_interf_cov_mats[kFreqBins][kNumMics][kNumMics];

    OsFloat wave_numbers[kFreqBins];
    OsFloat mask_thresholds[kFreqBins];
    OsFloat rxiws[kFreqBins];
    OsFloat rpsiws[kFreqBins];
    OsFloat reflected_rpsiws[kFreqBins];
    OsFloat new_mask[kFreqBins];
    OsFloat time_smooth_mask[kFreqBins];
    OsFloat final_mask[kFreqBins+1];

    // remove 50Hz and DC
    OsFloat coefs[6];
    OsFloat x_state[2];
    OsFloat y_state[2];
};

void PhaseAlignmentMasks(OsInt32 freq_bin,OsInt32 fft_size,OsInt32 sample_rate,
    OsFloat sound_speed,SensorPoint *inGeometry,OsFloat angle,complex *mat)
{
    // calculation steering vector
    OsFloat freq_in_hertz = 1.0 * freq_bin / fft_size * sample_rate;
    OsInt32 i;

    for(i = 0; i < kNumMics; ++i)
    {
        OsFloat distance = cos(angle)*inGeometry[i].x + sin(angle)*inGeometry[i].y;
        OsFloat phase_shift = -2 * M_PI * freq_in_hertz * distance / sound_speed;

        mat[i].re = cos(phase_shift);
        mat[i].im = sin(phase_shift);
    }
}

// 向量的内积
complex ConjugateDotProduct(complex *lhs,complex *rhs)
{
    complex result = init_c(0,0);
    OsInt32 i;

    for(i = 0; i < kNumMics; ++i)
    {
        complex temp = mul_c(rhs[i],conj_c(lhs[i]));
        result = add_c(result,temp);
    }
    return result;
}

void InitDelaySumMasks(Beamformer *inBeamformer)
{
    OsInt32 i,j;

    for (i = 0; i < kFreqBins; ++i)
    {
        complex *delay_sum_masks = inBeamformer->delay_sum_masks[i];
        complex norm_factor = init_c(1,0);
        complex sum = init_c(0,0);
        OsFloat sum_abs_norm_factor = 0;

        PhaseAlignmentMasks(i,kFftSize,inBeamformer->fs,kSpeedOfSoundMeterSeconds,
            inBeamformer->geometry,inBeamformer->targetAngle,delay_sum_masks);

        sum = ConjugateDotProduct(delay_sum_masks,delay_sum_masks);
        sum = sqrt_c(sum,2);
        norm_factor = div_c(norm_factor,sum);
        
        for(j = 0; j < kNumMics; ++j)
        {
            delay_sum_masks[j] = mul_c(delay_sum_masks[j],norm_factor);
            sum_abs_norm_factor += abs_c(delay_sum_masks[j]);
        }
        sum_abs_norm_factor = 1.0/sum_abs_norm_factor;
        for(j = 0; j < kNumMics; ++j)
        {
            inBeamformer->normalized_delay_sum_masks[i][j] = scale_c(delay_sum_masks[j],sum_abs_norm_factor);
        }
    }
}

void InitTargetCovMats(Beamformer *inBeamformer)
{
    complex numerator = init_c(1,0);
    OsInt32 i,j,k;

    for(k = 0; k < kFreqBins; ++k)
    {
        complex *in = inBeamformer->delay_sum_masks[k];
        complex (*out)[kNumMics] = inBeamformer->target_cov_mats[k];
        complex normalization_factor = init_c(0,0);

        // 计算每个频点的转向向量协方差矩阵
        for(i = 0; i < kNumMics; ++i)
        {
            for(j = 0; j < kNumMics; ++j)
            {
                out[i][j] = mul_c(in[i],conj_c(in[j]));
            }
        }
        // 取其迹
        for(i = 0; i < kNumMics; ++i)
        {
            normalization_factor = add_c(normalization_factor,out[i][i]);
        }
        normalization_factor = div_c(numerator,normalization_factor);
        for (i = 0; i < kNumMics; ++i)
        {
            for(j = 0; j < kNumMics; ++j)
            {
                out[i][j] = mul_c(out[i][j],normalization_factor);
            }
        }
    }
}

void UniformCovarianceMatrix(OsFloat wave_number,SensorPoint *inGeometry,complex (*mat)[kNumMics])
{
    OsInt32 i,j;

    for(i = 0; i < kNumMics; ++i)
    {
        for(j = 0; j < kNumMics; ++j)
        {
            OsFloat v;
            if(wave_number > 0.f)
            {
                // wave_number * distance得到的是频域的角频率
                v = j0(wave_number * Distance(inGeometry[i],inGeometry[j]));
                mat[i][j] = init_c(v,0);
            }
            else
            {
                v = (i == j) ? 1.0f : 0.0f;
                mat[i][j] = init_c(v,0);
            }
        }
    }
}

void AngledCovarianceMatrix(OsInt32 freq_bin,OsInt32 fs,OsFloat angle,SensorPoint *inGeometry,complex (*mat)[kNumMics])
{
    complex interf_cov_vector[kNumMics];
    OsInt32 i,j;

    PhaseAlignmentMasks(freq_bin,kFftSize,fs,kSpeedOfSoundMeterSeconds,inGeometry,angle,interf_cov_vector);
    for(i = 0; i < kNumMics; ++i)
    {
        complex x = interf_cov_vector[i];
        for(j = 0; j < kNumMics; ++j)
        {
            mat[i][j] = mul_c(x,conj_c(interf_cov_vector[j]));
        }
    }
}

void InitInterfConvMats(Beamformer *inBeamformer)
{
    complex uniform_cov_mat[kNumMics][kNumMics] = {0};
    complex angled_cov_mat[kNumMics][kNumMics] = {0};
    complex numerator = init_c(1,0);
    complex avg1,avg2;
    OsInt32 i,j,k;

    for(k = 0; k < kFreqBins; ++k)
    {
        complex normalization_factor = init_c(0,0);
        // 均匀干扰协方差矩阵uniform_cov_mat
        UniformCovarianceMatrix(inBeamformer->wave_numbers[k],inBeamformer->geometry,uniform_cov_mat);
        // 点源干扰协方差矩阵angled_cov_mat
        AngledCovarianceMatrix(k,inBeamformer->fs,inBeamformer->interfAngle,inBeamformer->geometry,angled_cov_mat);

        // Normalize matrices before averaging them
        for(i = 0; i < kNumMics; ++i)
        {
            normalization_factor = add_c(normalization_factor,uniform_cov_mat[i][i]);
        }
        normalization_factor = div_c(numerator,normalization_factor);
        for(i = 0; i < kNumMics; ++i)
        for(j = 0; j < kNumMics; ++j)
        {
            uniform_cov_mat[i][j] = mul_c(uniform_cov_mat[i][j],normalization_factor);
        }
        normalization_factor = init_c(0,0);
        for(i = 0; i < kNumMics; ++i)
        {
            normalization_factor = add_c(normalization_factor,angled_cov_mat[i][i]);
        }
        normalization_factor = div_c(numerator,normalization_factor);
        for(i = 0; i < kNumMics; ++i)
        for(j = 0; j < kNumMics; ++j)
        {
            angled_cov_mat[i][j] = mul_c(angled_cov_mat[i][j],normalization_factor);
        }
        // Average matrices
        avg1 = init_c(1-kBalance,0);
        avg2 = init_c(kBalance,0);
        for(i = 0; i < kNumMics; ++i)
        for(j = 0; j < kNumMics; ++j)
        {
            // 联合均匀干扰和点源干扰
            uniform_cov_mat[i][j] = mul_c(uniform_cov_mat[i][j],avg1);
            angled_cov_mat[i][j]  = mul_c(angled_cov_mat[i][j],avg2);

            inBeamformer->interf_cov_mats[k][i][j] = add_c(uniform_cov_mat[i][j],angled_cov_mat[i][j]);
            inBeamformer->reflected_interf_cov_mats[k][i][j] = conj_c(inBeamformer->interf_cov_mats[k][i][j]);
        }
    }
}

// Does conjugate(|norm_mat|) * |mat| * transpose(|norm_mat|),and norm_mat is a row vector
OsFloat Norm(complex (*mat)[kNumMics],complex *norm_mat)
{
    complex first_product = init_c(0,0);
    complex second_product = init_c(0,0);
    OsInt32 i,j;

    for (i = 0; i < kNumMics; ++i)
    {
        for (j = 0; j < kNumMics; ++j)
        {
            complex x = mul_c(conj_c(norm_mat[j]),mat[j][i]);
            first_product = add_c(first_product,x);
        }
        second_product = add_c(second_product,mul_c(first_product,norm_mat[i]));
        first_product = init_c(0,0);
    }
    return max(second_product.re,0.0f);
}

Beamformer* IcmBeamformerOpen(OsFloat (*inGeometry)[3],OsInt32 inFs,OsInt32 inFrameSize)
{
    OsFloat min_mic_spacing,kAliasingFreqHz,kHighMeanStartHz,kHighMeanEndHz;
    OsFloat targetAngle = M_PI/2;
    OsFloat interfAngle = M_PI/2-targetAngle;    // 120*M_PI/180.0;
    OsInt32 i;

    Beamformer *pBeamformer = (Beamformer*)icm_alloc(1,sizeof(Beamformer));
    assert(0 != pBeamformer);
    // 阵列形状赋值
    CenteredArray(inGeometry,kNumMics);
    for(i = 0; i < kNumMics; ++i)
    {
        pBeamformer->geometry[i].x = inGeometry[i][1];
        pBeamformer->geometry[i].y = inGeometry[i][0];
        pBeamformer->geometry[i].z = inGeometry[i][2];
    }
    min_mic_spacing     = GetMinimumSpacing(pBeamformer->geometry,kNumMics);
    kAliasingFreqHz     = kSpeedOfSoundMeterSeconds/(min_mic_spacing * (1.0 + abs(cos(targetAngle))));
    kHighMeanStartHz    = min(0.5f * kAliasingFreqHz,inFs / 2.f);
    kHighMeanEndHz      = min(0.75f * kAliasingFreqHz,inFs / 2.f);

    pBeamformer->fs                     = inFs;
    pBeamformer->frmsize                = inFrameSize;
    pBeamformer->targetAngle            = targetAngle;
    pBeamformer->interfAngle            = interfAngle;
    pBeamformer->low_mean_start_bin     = round(1.0 * kLowMeanStartHz * kFftSize / inFs);
    pBeamformer->low_mean_end_bin       = round(1.0 * kLowMeanEndHz * kFftSize / inFs);
    pBeamformer->high_mean_start_bin    = round(kHighMeanStartHz * kFftSize / inFs);
    pBeamformer->high_mean_end_bin      = round(kHighMeanEndHz * kFftSize / inFs);
    pBeamformer->num_output_channels    = 1;
    pBeamformer->init_delay             = kFftSize - gcd(pBeamformer->frmsize,kFftSize/2);
    pBeamformer->frame_offset           = 0;
    pBeamformer->shift_amount           = kFftSize/2;

    for(i = 0; i < kFreqBins; ++i)
    {
        pBeamformer->time_smooth_mask[i]    = 1.0;
        pBeamformer->final_mask[i]          = 1.0;
        pBeamformer->wave_numbers[i]        = 2.0 * M_PI * i / kFftSize * inFs / kSpeedOfSoundMeterSeconds;
        pBeamformer->mask_thresholds[i]     = kNumMics * kNumMics * kBeamwidthConstant * pBeamformer->wave_numbers[i] * pBeamformer->wave_numbers[i];
    }

    pBeamformer->window = (OsFloat*)icm_alloc(kFftSize,sizeof(OsFloat));
    assert(0 != pBeamformer->window);
    kaiser_c(kFftSize,2.0,pBeamformer->window);
    for(i = 0; i < kNumMics; ++i)
    {
        pBeamformer->x_buf[i]       = (OsFloat*)icm_alloc(kFftSize,sizeof(OsFloat));
        pBeamformer->in_block[i]    = (complex*)icm_alloc(kFreqBins,sizeof(complex));
        pBeamformer->norm_in_buf[i] = (OsFloat*)icm_alloc(pBeamformer->frmsize,sizeof(OsFloat));
        pBeamformer->in_buf[i]      = WebRtc_CreateBuffer(pBeamformer->frmsize+pBeamformer->init_delay,sizeof(OsFloat));
        // webrtc中使用的是init_delay，这里为了方便DSP，改成了shift_amount
        WebRtc_MoveReadPtr(pBeamformer->in_buf[i],-pBeamformer->shift_amount);
    }
    pBeamformer->y_buf          = (OsFloat*)icm_alloc(kFftSize,sizeof(OsFloat));
    pBeamformer->norm_out_buf   = (OsFloat*)icm_alloc(pBeamformer->frmsize,sizeof(OsFloat));
    pBeamformer->out_buf        = (OsFloat*)icm_alloc(pBeamformer->frmsize+pBeamformer->init_delay,sizeof(OsFloat));
    pBeamformer->out_block      = (complex*)icm_alloc(kFreqBins,sizeof(complex));

    pBeamformer->fft = icm_fft_init(kFftSize);
    assert(0 != pBeamformer->fft);

    InitDelaySumMasks(pBeamformer);
    InitTargetCovMats(pBeamformer);
    InitInterfConvMats(pBeamformer);

    for(i = 0; i < kFreqBins; ++i)
    {
        pBeamformer->rxiws[i]               = Norm(pBeamformer->target_cov_mats[i],pBeamformer->delay_sum_masks[i]);
        pBeamformer->rpsiws[i]              = Norm(pBeamformer->interf_cov_mats[i],pBeamformer->delay_sum_masks[i]);
        pBeamformer->reflected_rpsiws[i]    = Norm(pBeamformer->reflected_interf_cov_mats[i],pBeamformer->delay_sum_masks[i]);
    }

    {
        OsFloat K = tan(M_PI * 55.0 / inFs);
        OsFloat K2 = K * K;
        OsFloat Q = 0.75;
        OsFloat factor = 1.0 / (Q + K + K2 * Q);

        pBeamformer->coefs[0] = Q * factor;
        pBeamformer->coefs[1] = -2 * Q * factor;
        pBeamformer->coefs[2] = Q * factor;
        pBeamformer->coefs[3] = 1.0;
        pBeamformer->coefs[4] = 2 * Q * (K2 - 1) * factor;
        pBeamformer->coefs[5] = (K2 * Q - K + Q) * factor;
    }
    return pBeamformer;
}

void IcmBeamformerSetParamer(Beamformer *inBeamformer,OsFloat inTargetAngle,OsFloat inInterfAngle)
{
    OsFloat min_mic_spacing     = GetMinimumSpacing(inBeamformer->geometry,kNumMics);
    OsFloat kAliasingFreqHz     = kSpeedOfSoundMeterSeconds/(min_mic_spacing * (1.0 + abs(cos(inTargetAngle))));
    OsFloat kHighMeanStartHz    = min(0.5f * kAliasingFreqHz,inBeamformer->fs/2.f);
    OsFloat kHighMeanEndHz      = min(0.75f * kAliasingFreqHz,inBeamformer->fs/2.f);
    OsInt32 i;

    inBeamformer->targetAngle   = inTargetAngle;
    inBeamformer->interfAngle   = inInterfAngle;

    InitDelaySumMasks(inBeamformer);
    InitTargetCovMats(inBeamformer);
    InitInterfConvMats(inBeamformer);

    for(i = 0; i < kFreqBins; ++i)
    {
        inBeamformer->rxiws[i]              = Norm(inBeamformer->target_cov_mats[i],inBeamformer->delay_sum_masks[i]);
        inBeamformer->rpsiws[i]             = Norm(inBeamformer->interf_cov_mats[i],inBeamformer->delay_sum_masks[i]);
        inBeamformer->reflected_rpsiws[i]   = Norm(inBeamformer->reflected_interf_cov_mats[i],inBeamformer->delay_sum_masks[i]);
    }
}

void IcmBeamformerClose(Beamformer *inBeamformer)
{
    OsInt32 i;

    for (i = 0; i < kNumMics; ++i)
    {
        icm_free(inBeamformer->norm_in_buf[i]);
        WebRtc_FreeBuffer(inBeamformer->in_buf[i]);
        icm_free(inBeamformer->x_buf[i]);
        icm_free(inBeamformer->in_block[i]);
    }
    icm_free(inBeamformer->y_buf);
    icm_free(inBeamformer->norm_out_buf);
    icm_free(inBeamformer->out_buf);
    icm_free(inBeamformer->out_block);
    icm_fft_uninit(inBeamformer->fft);
    icm_free(inBeamformer->window);
    icm_free(inBeamformer);
}

OsFloat MaskRangeMean(OsFloat *time_smooth_mask,OsInt32 first,OsInt32 last)
{
    OsFloat sum = 0;
    OsInt32 i;
    
    for(i = first; i < last; ++i)
    {
        sum += time_smooth_mask[i];
    }
    return sum/(last-first);
}

OsFloat CalculatePostfilterMask(complex (*cov_mat)[kNumMics],complex *eig_m,
    OsFloat rpsiw,OsFloat ratio_rxiw_rxim,OsFloat rmw_r,OsFloat mask_threshold)
{
    // 这里的cov_mat是interf_cov_mats或者reflected_interf_cov_mats
    OsFloat rpsim = Norm(cov_mat,eig_m);
    OsFloat ratio = 0;
    OsFloat numerator;
    OsFloat denominator;
    OsFloat mask = 1.0f;

    if (rpsim > 0.f)
    {
        // rpsiw是rpsiws或者reflected_rpsiws
        ratio = rpsiw/rpsim;
    }
    // rmw_r是eig_m与delay_sum_masks共轭的内积的功率
    numerator = rmw_r - ratio;
    // ratio_rxiw_rxim是：rxiws除以eig_m与target_cov_mats的范数rxim
    denominator = ratio_rxiw_rxim - ratio;

    if (denominator > mask_threshold)
    {
        OsFloat lambda = numerator / denominator;
        mask = max(lambda * ratio_rxiw_rxim / rmw_r,kMaskMinimum);
    }
    return mask;
}

void ProcessAudio(Beamformer *inBeamformer,complex **input,complex *output)
{
    OsFloat *new_mask           = inBeamformer->new_mask;
    OsFloat *final_mask         = inBeamformer->final_mask;
    OsFloat *time_smooth_mask   = inBeamformer->time_smooth_mask;
    OsInt32 low_mean_start_bin  = inBeamformer->low_mean_start_bin;
    OsInt32 low_mean_end_bin    = inBeamformer->low_mean_end_bin;
    OsInt32 high_mean_start_bin = inBeamformer->high_mean_start_bin;
    OsInt32 high_mean_end_bin   = inBeamformer->high_mean_end_bin;
    OsFloat low_frequency_mask,high_pass_postfilter_mask;
    complex eig_m[kNumMics];
    OsInt32 i,j;

    for (i = low_mean_start_bin; i <= high_mean_end_bin; ++i)
    {
        OsFloat eig_m_norm_factor = 0;
        OsFloat ratio_rxiw_rxim = 0;
        OsFloat rxim;
        OsFloat sum_squares = 0;
        OsFloat abs_value;
        OsFloat rmw,rmw_r;

        // 输入信号归一化 
        for (j = 0; j < kNumMics; ++j)
        {
            eig_m[j] = input[j][i];

            abs_value = abs_c(eig_m[j]);
            sum_squares += abs_value*abs_value;
        }
        eig_m_norm_factor = sqrt(sum_squares);
        if (sum_squares != 0.0f)
        {
            eig_m_norm_factor = 1.0/eig_m_norm_factor;
            for(j = 0; j < kNumMics; ++j)
            {
                eig_m[j] = scale_c(eig_m[j],eig_m_norm_factor);
            }
        }

        rxim = Norm(inBeamformer->target_cov_mats[i],eig_m);
        if (rxim > 0.f)
        {
            ratio_rxiw_rxim = inBeamformer->rxiws[i]/rxim;
        }

        rmw = abs_c(ConjugateDotProduct(inBeamformer->delay_sum_masks[i],eig_m));
        rmw_r = rmw*rmw;
        
        new_mask[i]  = CalculatePostfilterMask(inBeamformer->interf_cov_mats[i],eig_m,
            inBeamformer->rpsiws[i],ratio_rxiw_rxim,rmw_r,inBeamformer->mask_thresholds[i]);
        new_mask[i] *= CalculatePostfilterMask(inBeamformer->reflected_interf_cov_mats[i],eig_m,
            inBeamformer->reflected_rpsiws[i],ratio_rxiw_rxim,rmw_r,inBeamformer->mask_thresholds[i]);
    }

    for (i = low_mean_start_bin; i <= high_mean_end_bin; ++i)
    {
        time_smooth_mask[i] = kMaskTimeSmoothAlpha * new_mask[i] + (1 - kMaskTimeSmoothAlpha) * time_smooth_mask[i];
    }

    low_frequency_mask = MaskRangeMean(time_smooth_mask,low_mean_start_bin,low_mean_end_bin+1);
    for(i = 0; i < low_mean_start_bin; ++i)
    {
        time_smooth_mask[i] = low_frequency_mask;
    }
    high_pass_postfilter_mask = MaskRangeMean(time_smooth_mask,high_mean_start_bin,high_mean_end_bin+1);
    for(i = high_mean_end_bin+1; i < kFreqBins; ++i)
    {
        time_smooth_mask[i] = high_pass_postfilter_mask;
    }
    for(i = 0; i < kFreqBins; ++i)
    {
        final_mask[i] = time_smooth_mask[i];
    }
    final_mask[kFreqBins] = final_mask[kFreqBins-1];    // 多出一个，防止因平滑操作导致数值不稳定
    for (i = low_mean_start_bin; i < kFreqBins; ++i)
    {
        final_mask[i] = kMaskFrequencySmoothAlpha * final_mask[i] + (1-kMaskFrequencySmoothAlpha) * final_mask[i-1];
    }
    for (i = high_mean_end_bin; i >= 0; --i)
    {
        final_mask[i] = kMaskFrequencySmoothAlpha * final_mask[i] + (1-kMaskFrequencySmoothAlpha) * final_mask[i+1];
    }

    for(i = 0; i < kFreqBins; ++i)
    {
        complex *delay_sum_mask = inBeamformer->normalized_delay_sum_masks[i];
        output[i] = init_c(0,0);

        for(j = 0; j < kNumMics; ++j)
        {
            output[i] = add_c(output[i],mul_c(input[j][i],delay_sum_mask[j]));
        }
output[i] = scale_c(output[i],final_mask[i]);
    }
}

void ProcessBlock(Beamformer *inBeamformer,OsFloat **input,OsFloat *output)
{
    OsInt32 frmsize     = inBeamformer->frmsize;
    OsFloat *window     = inBeamformer->window;
    OsFloat *out_buf    = inBeamformer->out_buf;
    OsFloat *y_buf      = inBeamformer->y_buf;
    OsInt32 first_frame_in_block;
    OsInt32 i,j;

    for(i = 0; i < kNumMics; ++i)
    {
        WebRtc_WriteBuffer(inBeamformer->in_buf[i],input[i],frmsize);
    }
    first_frame_in_block = inBeamformer->frame_offset;
    
    while(first_frame_in_block < frmsize)
    {
        for (i = 0; i < kNumMics; ++i)
        {
            WebRtc_ReadBuffer(inBeamformer->in_buf[i],NULL,inBeamformer->x_buf[i],kFftSize);
            WebRtc_MoveReadPtr(inBeamformer->in_buf[i],-kFftSize+inBeamformer->shift_amount);
            for(j = 0; j < kFftSize; ++j)
            {
                inBeamformer->x_buf[i][j] = inBeamformer->x_buf[i][j] * window[j];
            }
            icm_fft(inBeamformer->fft,inBeamformer->x_buf[i],inBeamformer->in_block[i]);
        }

        ProcessAudio(inBeamformer,inBeamformer->in_block,inBeamformer->out_block);

        icm_ifft(inBeamformer->fft,inBeamformer->out_block,y_buf);
        for(j = 0; j < kFftSize; ++j)
        {
            y_buf[j] = y_buf[j] * window[j];
        }
        for(j = 0; j < kFftSize; ++j)
        {
            out_buf[j+first_frame_in_block] += y_buf[j];
        }
        first_frame_in_block += inBeamformer->shift_amount;
    }
    memcpy(output,out_buf,frmsize*sizeof(OsFloat));
    memcpy(out_buf,out_buf+frmsize,inBeamformer->init_delay*sizeof(OsFloat));
    memset(out_buf+inBeamformer->init_delay,0,frmsize*sizeof(OsFloat));
    inBeamformer->frame_offset = first_frame_in_block - frmsize;
}

void IcmBeamformerProcess(Beamformer *inBeamformer,OsInt16 *inPcm,OsInt32 inLen,OsInt16 *ioPcm,OsInt32 *ioLen)
{
    OsFloat *norm_out_buf = inBeamformer->norm_out_buf;
    OsInt32 frmsize = inBeamformer->frmsize;
    OsFloat *norm_in_buf[kNumMics];
    OsInt32 i,j;

    // PCM数据归一化
    for(i = 0; i < kNumMics; ++i)
    {
        OsInt32 idx = i;
        norm_in_buf[i] = inBeamformer->norm_in_buf[i];
        for(j = 0; j < frmsize; ++j)
        {
            norm_in_buf[i][j] = inPcm[idx]/S16_MAX;
            idx += kNumMics;
        }
    }

    ProcessBlock(inBeamformer,norm_in_buf,norm_out_buf);
    //{
    //    // remove 50Hz and DC
    //    OsFloat *x = &inBeamformer->x_state[1];
    //    OsFloat *y = &inBeamformer->y_state[1];
    //    OsFloat *b = inBeamformer->coefs;
    //    OsFloat *a = b+3;
    //    OsInt32 i;

    //    for(i = 0; i < frmsize; i++)
    //    {
    //        OsFloat x0 = norm_out_buf[i];
    //        OsFloat y_new = b[0]*x0 + b[1]*x[0] + b[2]*x[-1] - a[1]*y[0] - a[2]*y[-1];

    //        x[-1] = x[0]; x[0] = x0;
    //        y[-1] = y[0]; y[0] = y_new;
    //        x0 = y_new;

    //        norm_out_buf[i] = y_new;
    //    }
    //}

    // 归一化数据转为PCM数据
    for(i = 0; i < frmsize; ++i)
    {
        OsFloat v = sat(norm_out_buf[i],-1,1);
        ioPcm[i] = round(v*S16_MAX);
    }



    *ioLen = frmsize;
}