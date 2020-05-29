#include "sbc.h"
#include <android/log.h>
#include <jni.h>

#define  LOG_TAG    "libsbc"
#define  LOGI(...)  __android_log_print(ANDROID_LOG_INFO, LOG_TAG, __VA_ARGS__)
//#define PCM_SIZE1 512
//#define PCM_SIZE4 (PCM_SIZE1 * 4)
//#define SBC_SIZE1 123
//#define SBC_SIZE4 (SBC_SIZE1 * 4)
//#define packetHeadLen 6

sbc_t sbc;
sbc_t sbc_encode_t;
int codesize;
int codesize_encode;

JNIEXPORT void JNICALL Java_com_sabinetek_swiss_jni_codec_CodecSbc_SbcInit(JNIEnv *env,
                                                                           jobject cls,
                                                                           jint flags) {

    LOGI("SbcInit: %d", flags);
    sbc_init(&sbc, flags);
    codesize = sbc_get_codesize(&sbc);

    sbc_init(&sbc_encode_t, flags);
    codesize_encode = sbc_get_codesize(&sbc_encode_t);

    sbc_encode_t.frequency = SBC_FREQ_48000;
    sbc_encode_t.subbands = SBC_SB_8;
    sbc_encode_t.mode = SBC_MODE_JOINT_STEREO;
    sbc_encode_t.endian = SBC_LE;
    sbc_encode_t.bitpool = 40;//55;
    sbc_encode_t.allocation = SBC_AM_LOUDNESS;
    sbc_encode_t.blocks = SBC_BLK_16;


}

JNIEXPORT jint JNICALL Java_com_sabinetek_swiss_jni_codec_CodecSbc_SbcDecode(JNIEnv *env,
                                                                             jclass cls,
                                                                             jbyteArray in,
                                                                             jint inlen,
                                                                             jbyteArray out) {
    jbyte *jinput = (*env)->GetByteArrayElements(env, in, JNI_FALSE);
    jbyte *joutput = (*env)->GetByteArrayElements(env, out, JNI_FALSE);

    int output_len = (*env)->GetArrayLength(env, out);
    unsigned char *inp = (unsigned char *) jinput;
    unsigned char *outp = (unsigned char *) joutput;
    size_t decoded;


    int framelen = 0;
    int decoded_pcm_len = 0;
    int output_buf_size = 0;
    int input_buf_size = 0;
    unsigned int frame_decoded_pcm_len = 0;
    unsigned int decoded_sbc_len = 0;

    output_buf_size = output_len;
    input_buf_size = inlen;

    framelen = sbc_decode(&sbc, inp, input_buf_size, outp, output_buf_size, &frame_decoded_pcm_len);
    decoded_sbc_len += framelen;
    decoded_pcm_len += frame_decoded_pcm_len;
//    LOGI("framelen : %d", framelen);
//    LOGI("decoded_pcm_len : %d", frame_decoded_pcm_len);

    while (framelen > 0) {
        if (decoded_pcm_len > output_buf_size)
            break;

        /* push the pointer in the file forward to the next bit to be
         * decoded tell the decoder to decode up to the remaining
         * length of the file (!) */


        framelen = sbc_decode(&sbc, (char *) inp + decoded_sbc_len,
                              input_buf_size - decoded_sbc_len,
                              (char *) outp + decoded_pcm_len, output_buf_size - decoded_pcm_len,
                              &frame_decoded_pcm_len);
//        LOGI("while framelen : %d", framelen);
//        LOGI("while decoded_pcm_len : %d", frame_decoded_pcm_len);
        if (framelen < 0)
            break;
        decoded_sbc_len += framelen;
        decoded_pcm_len += frame_decoded_pcm_len;


    }

    (*env)->ReleaseByteArrayElements(env, in, jinput, JNI_FALSE);
    (*env)->ReleaseByteArrayElements(env, out, joutput, JNI_FALSE);
//    LOGI("decoded_pcm_len: %d", decoded_pcm_len);
    return decoded_pcm_len;
}

JNIEXPORT jint JNICALL Java_com_sabinetek_swiss_jni_codec_CodecSbc_SbcEncode(JNIEnv *env,
                                                                             jclass cls,
                                                                             jbyteArray in,
                                                                             jint inlen,
                                                                             jbyteArray out) {
    jbyte *jinput = (*env)->GetByteArrayElements(env, in, JNI_FALSE);
    jbyte *joutput = (*env)->GetByteArrayElements(env, out, JNI_FALSE);

    int outlen = (*env)->GetArrayLength(env, out);

    unsigned char *inp = (unsigned char *) jinput;
    unsigned char *outp = (unsigned char *) joutput;

    int nframes = 0;
    int len_res = 0;
    int encoded_pcm_len = 0;
    int output_buf_size = 0;
    int input_buf_size = 0;
    unsigned int frame_encoded_pcm_len = 0;
    int frame_encoded_sbc_len = 0;
    unsigned int encoded_sbc_len = 0;

    output_buf_size = outlen;
    input_buf_size = inlen;


    frame_encoded_pcm_len = sbc_get_codesize(&sbc_encode_t);

    //according to the framesize;
    if (input_buf_size % frame_encoded_pcm_len != 0)
        return -1;

    nframes = input_buf_size / frame_encoded_pcm_len;

    encoded_pcm_len = 0;
    encoded_sbc_len = 0;


    while (encoded_pcm_len < input_buf_size) {
        len_res = sbc_encode(&sbc_encode_t, (char *) inp + encoded_pcm_len, frame_encoded_pcm_len,
                             (char *) outp + encoded_sbc_len, output_buf_size - encoded_sbc_len,
                             &frame_encoded_sbc_len);
        if (len_res != frame_encoded_pcm_len || frame_encoded_sbc_len <= 0) {
            return -2;
        }
        encoded_pcm_len += len_res;
        encoded_sbc_len += frame_encoded_sbc_len;
    }

    (*env)->ReleaseByteArrayElements(env, in, jinput, JNI_FALSE);
    (*env)->ReleaseByteArrayElements(env, out, joutput, JNI_FALSE);

    return encoded_sbc_len;
}

JNIEXPORT void JNICALL Java_com_sabinetek_swiss_jni_codec_CodecSbc_SbcExit(JNIEnv *env,
                                                                           jclass cls) {
    sbc_finish(&sbc);
}
