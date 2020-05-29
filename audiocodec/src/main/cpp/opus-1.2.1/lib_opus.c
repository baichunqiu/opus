#include <string.h>
#include <stdio.h>
#include <jni.h>
#include <android/log.h>
#include <opus.h>

#define TAG "Opus_JNI"
#define LOGV(...) __android_log_print(ANDROID_LOG_VERBOSE, TAG,__VA_ARGS__)
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG , TAG,__VA_ARGS__)
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO  , TAG,__VA_ARGS__)
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN  , TAG,__VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR , TAG,__VA_ARGS__)
#define MAX_PACKET 1500

OpusEncoder *enc = NULL;
OpusDecoder *dec = NULL;

JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *vm, void *reserved) {
    LOGD(TAG, "load library");
    JNIEnv *env = NULL;
    jint result = -1;

    if ((*vm)->GetEnv(vm, (void **) &env, JNI_VERSION_1_4) != JNI_OK) {
        LOGE("tag", "load library error ");
        return JNI_ERR;
    }
    result = JNI_VERSION_1_4;
    LOGD(TAG, "load library result");
    return result;
}

JNIEXPORT void JNICALL JNI_OnUnload(JavaVM *vm, void *reserved) {
    __android_log_print(ANDROID_LOG_ERROR, "tag", "library was unload");
}

jint Java_com_sabinetek_swiss_jni_codec_CodecOpus_initEncode(JNIEnv *env, jobject thiz,
                                                             jint sampling_rate, jint channels,
                                                             jint bitrate, jint complexity,
                                                             jint bitsPerSample) {

    int err;
    int use_vbr;
    int cvbr = 0;
    int application = OPUS_APPLICATION_AUDIO;
    opus_int32 bitrate_bps = 0;
    int max_payload_bytes;
    int forcechannels;
    int variable_duration = OPUS_FRAMESIZE_ARG;

    /* defaults: */
    use_vbr = 0;
    max_payload_bytes = MAX_PACKET;
    forcechannels = OPUS_AUTO;

    if (enc != NULL) {
        LOGE("CELT encoder already initialed!!!");
        return -1;
    }
    sampling_rate = 48000;

    if ((sampling_rate != 48000) || (channels != 2) || (bitsPerSample != 16)) {
        LOGE("CELT encoder parameters wrong!!!");
        return -1;
    }

    enc = opus_encoder_create(sampling_rate, channels, application, &err);
    if (err != OPUS_OK) {
        LOGE("CELT encoder open fail!!!");
        return -1;
    }
    opus_encoder_ctl(enc, OPUS_SET_BITRATE(bitrate_bps));
    opus_encoder_ctl(enc, OPUS_SET_BANDWIDTH(OPUS_BANDWIDTH_FULLBAND));
    opus_encoder_ctl(enc, OPUS_SET_VBR(use_vbr));
    opus_encoder_ctl(enc, OPUS_SET_VBR_CONSTRAINT(cvbr));
    opus_encoder_ctl(enc, OPUS_SET_COMPLEXITY(complexity));
    opus_encoder_ctl(enc, OPUS_SET_FORCE_CHANNELS(forcechannels));
    opus_encoder_ctl(enc, OPUS_SET_LSB_DEPTH(16));
    opus_encoder_ctl(enc, OPUS_SET_EXPERT_FRAME_DURATION(variable_duration));
    opus_encoder_ctl(enc, OPUS_SET_SIGNAL(OPUS_SIGNAL_MUSIC));

    return 0;

}

jint Java_com_sabinetek_swiss_jni_codec_CodecOpus_initDecode(JNIEnv *env, jobject thiz,
                                                             jint sampling_rate, jint channels,
                                                             jint version) {

    int err;
    sampling_rate = 48000;

    if (dec != NULL) {
        LOGE("CELT decoder already initialed!!!");
        return -1;
    }
    if ((sampling_rate != 48000) || (channels != 2)) {
        LOGE("CELT decoder parameters wrong!!!");
        return -1;
    }

    dec = opus_decoder_create(sampling_rate, channels, &err);
    if (err != OPUS_OK) {
        LOGE("Cannot create decoder: %s\n", opus_strerror(err));

        return -1;
    }
    opus_decoder_ctl(dec, OPUS__SET_SABINE_VERSION(version));
    return 0;
}

static void inline int_to_char(opus_uint32 i, unsigned char ch[4]) {
    ch[0] = i >> 24;
    ch[1] = (i >> 16) & 0xFF;
    ch[2] = (i >> 8) & 0xFF;
    ch[3] = i & 0xFF;
}

static opus_uint32 inline char_to_int(unsigned char ch[4]) {
    return ((opus_uint32) ch[0] << 24) | ((opus_uint32) ch[1] << 16)
           | ((opus_uint32) ch[2] << 8) | (opus_uint32) ch[3];
}

jint Java_com_sabinetek_swiss_jni_codec_CodecOpus_encode(JNIEnv *env, jobject thisz,
                                                         jbyteArray in, jint frame_size,
                                                         jbyteArray out) {

    int len;
    unsigned char *pdata;
    opus_uint32 enc_final_range[2];
    unsigned char int_field[4];

    const opus_int16 *pcm;

    jbyte *pin = (*env)->GetByteArrayElements(env, in, JNI_FALSE);
    jbyte *pout = (*env)->GetByteArrayElements(env, out, JNI_FALSE);

    pcm = (opus_int16 *) pin;
    if (enc == NULL) {
        LOGE("CELT encoder NOT initialed!!!");
        return -1;
    }

    pdata = (unsigned char *) pout;
    pdata += 4;

    len = opus_encode(enc, pcm, frame_size / 4, pdata, 1500);

    opus_encoder_ctl(enc, OPUS_GET_FINAL_RANGE(&enc_final_range[0]));

    if (len < 0) {
        LOGE("CELT encoder error!!!");
        return -1;
    }

    pdata = (unsigned char *) pout;
    int_to_char(len, pdata);

    (*env)->ReleaseByteArrayElements(env, in, pin, JNI_FALSE);
    (*env)->ReleaseByteArrayElements(env, out, pout, JNI_FALSE);
    return len + 4;

}

jint Java_com_sabinetek_swiss_jni_codec_CodecOpus_decode(JNIEnv *env, jobject thisz,
                                                         jbyteArray dat, jint length,
                                                         jbyteArray pcm) {
    jint decodedSize;

    int i;
    int len;
    const unsigned char *data;
    unsigned char ch[4];
    opus_uint32 enc_final_range[2];
    opus_uint32 dec_final_range;
//	int output_samples = 2 * 48000;
    int pcm_frame_size_u16 = 2 * 960;
    int max_payload_bytes = MAX_PACKET;
    int decode_samples = 0, total_samples = 0;

    opus_int16 *out_ptr;

    jbyte *pdata = (*env)->GetByteArrayElements(env, dat, JNI_FALSE);
    jbyte *pout = (*env)->GetByteArrayElements(env, pcm, JNI_FALSE);

    out_ptr = (opus_int16 *) pout;
    len = length;

    if ((len > max_payload_bytes) || (len < 0)) {
        LOGE("Invalid payload length: %d\n", len);
    }

    data = (unsigned char *) pdata;
//	decode_samples = opus_decode(dec, data, len, out_ptr, output_samples, 0);

    opus_uint32 num, num1, num2, num3, num4;
    num1 = (opus_uint32) (*data) & 0x00 << 24;
    num2 = ((opus_uint32) *(data + 1)) & 0x00 << 16;
    num3 = ((opus_uint32) *(data + 2)) << 8;
    num4 = ((opus_uint32) *(data + 3));
    num = num1 + num2 + num3 + num4;
    int count = len / (num + 4);

    while (count > 0) {
        data += 4;
        decode_samples = opus_decode(dec, data, num, out_ptr,
                                     pcm_frame_size_u16, 0);
        data += num;
        out_ptr += decode_samples * 2;
        total_samples += decode_samples;
        count--;
    };

    (*env)->ReleaseByteArrayElements(env, dat, pdata, JNI_FALSE);
    (*env)->ReleaseByteArrayElements(env, pcm, pout, JNI_FALSE);
    return total_samples * 2 * 2;

}

void Java_com_sabinetek_swiss_jni_codec_CodecOpus_closeEncode(JNIEnv *env, jobject thiz) {

    if (enc == NULL)
        return;
    opus_encoder_destroy(enc);
    enc = NULL;
}

void Java_com_sabinetek_swiss_jni_codec_CodecOpus_closeDecode(JNIEnv *env, jobject thiz) {
    if (dec == NULL)
        return;
    dec = NULL;
}

