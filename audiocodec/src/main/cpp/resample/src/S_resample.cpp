#include "S_resample.h"
#include <android/log.h>
#define TAG    "S_resample"
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG,__VA_ARGS__)
S_resample::S_resample(uint32_t nNewSampleRate, uint32_t nOldSampleRate,
		uint32_t nBit, uint32_t nChannel, uint32_t nLengthSrc) {
	S_resample::m_kw_resample = new Mixer_RS_KW();
	S_resample::m_liner_resample = new Mixer_RS_Linear();
	S_resample::m_nBit = nBit;
	S_resample::m_nChannel = nChannel;
	S_resample::m_nNewSampleRate = nNewSampleRate;
	S_resample::m_nOldSampleRate = nOldSampleRate;
	S_resample::m_nLengthSrc = nLengthSrc;
	LOGI("init...%d -> %d", nOldSampleRate, nNewSampleRate);
	s_resample_ioctl();
}
bool S_resample::s_resample_reset() {
	if (!(m_kw_resample->Reset())) {
		return false;
	}
	if (!(m_liner_resample->Reset())) {
		return false;
	}
	return true;
}
bool S_resample::s_resample_ioctl() {
	if (m_nBit != 16) {
		return false;
	}
	if ((m_nChannel != 1) && (m_nChannel != 2)) {
		return false;
	}
	if (!(m_liner_resample->SetMethod(0, m_nChannel))) {
		return false;
	}
	if (!(m_liner_resample->SetRatio(m_nNewSampleRate, m_nOldSampleRate))) {
		return false;
	};
	if (!(m_kw_resample->SetMethod(0, m_nChannel))) {
		return false;
	}
	if (!(m_kw_resample->SetRatio(m_nNewSampleRate, m_nOldSampleRate))) {
		return false;
	};

//	LOGI("ioct1...");
	return true;

}

bool S_resample::s_resample_execute(LIBPCM_RESAMPLE_EXE_SRC_PARAM_T &src_param,
		LIBPCM_RESAMPLE_EXE_DES_PARAM_T &des_param) {

	int32_t src_adjust_len;
	uint32_t per_sample;
	int32_t tmp_src_len, tmp_src_used_len, tmp_des_len, tmp_des_used_len;
	if (((src_param.src_data == NULL) || (des_param.des_data == NULL))) {
		return false;
	}
	if ((src_param.src_data_len == 0) || (des_param.des_data_len == 0)) {
		src_param.src_used_data_len = 0;
		des_param.des_used_data_len = 0;
		LOGI("execute...0000");
		return true;
	}
	///< if new sample rate is equal with old sample rate

	if (m_nNewSampleRate == m_nOldSampleRate) {
		///< memcpy
		if (src_param.src_data_len < des_param.des_data_len) {
			memcpy(des_param.des_data, src_param.src_data,
					src_param.src_data_len);
			src_param.src_used_data_len = src_param.src_data_len;
			des_param.des_used_data_len = src_param.src_data_len;
		} else {
			memcpy(des_param.des_data, src_param.src_data,
					des_param.des_data_len);
			src_param.src_used_data_len = des_param.des_data_len;
			des_param.des_used_data_len = des_param.des_data_len;
		}
		return true;
	}
	///< calculate the src_data_len base on des_data_len, update src_data_len
	src_adjust_len = (int32_t) ((int64_t) des_param.des_data_len
			* (int64_t) m_nOldSampleRate / m_nNewSampleRate);
	//src_adjust_len = src_param.src_data_len;
	///< the output buffer is big enough, can hold all the data
	if (src_param.src_data_len < src_adjust_len) {
		src_adjust_len = src_param.src_data_len;
	}
	src_adjust_len = src_adjust_len & 0xfffffffc;
	if (src_adjust_len == 0) {
		src_param.src_used_data_len = 0;
		des_param.des_used_data_len = 0;
		return true;
	}
	if (m_nNewSampleRate > m_nOldSampleRate) {
		if (src_adjust_len > m_nLengthSrc) {
			BOOL ret;
			ret = S_resample::m_kw_resample->SetPCMBufLen(src_adjust_len + 256);
			if (!ret) {
				return false;
			}
			m_nLengthSrc = src_adjust_len;
		}
	}
	///< bytes to samples convert
	per_sample = m_nBit * m_nChannel / 8;
	tmp_src_len = src_adjust_len / per_sample;
	tmp_des_len = des_param.des_data_len / per_sample;
//	LOGI("execute...1111 %d %d %d %d", src_param.src_data_len, src_param.src_used_data_len, des_param.des_data_len, des_param.des_used_data_len);

	if (S_resample::m_nNewSampleRate > S_resample::m_nOldSampleRate) {
		//LOGI("execute...pre1%);
		//LOGI("execute...pre2%D"+tmp_src_len);

		tmp_des_used_len = m_kw_resample->ConvertUp(
				(short int*) des_param.des_data, tmp_des_len,
				(short int*) src_param.src_data, tmp_src_len,
				&tmp_src_used_len);
		/*tmp_des_used_len = m_liner_resample->ConvertUp((short int*)des_param.des_data, tmp_des_len, \
			(short int*)src_param.src_data, tmp_src_len, &tmp_src_used_len);*/

		if (tmp_des_used_len == -1) {
			return false;
		}
	} else {
		tmp_des_used_len = S_resample::m_liner_resample->ConvertDown(
				(short int*) des_param.des_data, tmp_des_len,
				(short int*) src_param.src_data, tmp_src_len,
				&tmp_src_used_len);
		if (tmp_des_used_len == -1) {
			return false;
		}
	}
	///< samples to bytes convert
	src_param.src_used_data_len = tmp_src_len * per_sample;
	des_param.des_used_data_len = tmp_des_used_len * per_sample;
	return true;
}
