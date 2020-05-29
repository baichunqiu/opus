//
// Created by ztimc on 2017/12/20.
//
#include <jni.h>
#include <beamformer.h>
#include <android/log.h>

#define TAG "log" // 这个是自定义的LOG的标识
#define LOGD(...) __android_log_print(ANDROID_LOG_DEBUG,TAG ,__VA_ARGS__) // 定义LOGD类型
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG ,__VA_ARGS__) // 定义LOGI类型
#define LOGW(...) __android_log_print(ANDROID_LOG_WARN,TAG ,__VA_ARGS__) // 定义LOGW类型
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR,TAG ,__VA_ARGS__) // 定义LOGE类型
#define LOGF(...) __android_log_print(ANDROID_LOG_FATAL,TAG ,__VA_ARGS__) // 定义LOGF类型
;

extern "C"
JNIEXPORT jlong JNICALL Java_com_sabinetek_swiss_jni_beamformer_Beanformer_open(JNIEnv *env,
                                                                              jobject obj,
                                                                              jint sampleRate,
                                                                              jint frameSize) {

    OsFloat mic_pos[2][3] = {{0, 0,    0},
                             {0, 0.17, 0}};

    Beamformer *pBeamformer = IcmBeamformerOpen(mic_pos, sampleRate, frameSize);

    return reinterpret_cast<jlong>(pBeamformer);
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_com_sabinetek_swiss_jni_beamformer_Beanformer_nativeProcess(JNIEnv *env,
                                                               jobject obj,
                                                               jlong nativePointer,
                                                               jbyteArray inPcm,
                                                               jint inLen) {
    jbyte *inData = (jbyte *) malloc((size_t) inLen);
    env->GetByteArrayRegion(inPcm, 0, inLen, inData);

    OsInt16 *outData = (OsInt16 *) malloc(5000);
    OsInt32 outSize;

    IcmBeamformerProcess(reinterpret_cast<Beamformer *>(nativePointer),
                         (OsInt16 *) inData,
                         inLen,
                         outData,
                         &outSize); 

    jbyteArray result = env->NewByteArray(outSize * 2);
    env->SetByteArrayRegion(result, 0, outSize * 2, (const jbyte *) outData);
    free(inData);
    free(outData);

    return result;
}


extern "C"
JNIEXPORT void JNICALL Java_com_sabinetek_swiss_jni_beamformer_Beanformer_nativeClose(JNIEnv *env,
                                                                                    jobject obj,
                                                                                    jlong nativePointer) {
    if (nativePointer == 0) return;
    Beamformer *beamformer = reinterpret_cast<Beamformer *>(nativePointer);
    IcmBeamformerClose(beamformer);
}