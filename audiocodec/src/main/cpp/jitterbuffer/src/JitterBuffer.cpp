#include <jni.h>
#include <stdio.h>
#include <android/log.h>
#include <DelayBuf.h>

#define TAG    "jitterbuffer"
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG,__VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR,TAG,__VA_ARGS__)


extern "C"
JNIEXPORT jlong JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_init(JNIEnv *env,
                                                  jobject jobj,
                                                  jint sampling_rate,
                                                  jint channel,
                                                  jint frame_time,
                                                  jint max_delay) {

    CDelayBuf *cDelayBuf = new CDelayBuf(sampling_rate, channel, frame_time, max_delay);

    return reinterpret_cast<jlong>(cDelayBuf);
}

extern "C"
JNIEXPORT jint JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_put(JNIEnv *env,
                                                 jobject jobj,
                                                 jlong native_pointer,
                                                 jbyteArray j_data,
                                                 jint j_length) {
    if (NULL == native_pointer) {
        LOGE("YOU MUST INIT BEFER USE");
        return -1;
    }
    CDelayBuf *delay_buffer = reinterpret_cast<CDelayBuf *>(native_pointer);
    jbyte *data = env->GetByteArrayElements(j_data,JNI_FALSE);
    int result = delay_buffer->DelayBufPut(reinterpret_cast<int16_t *>(data), j_length/2);
    if (result == -1) {
        LOGE("PUT ERROR");
        return -1;
    }
    return 0;

}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_get(JNIEnv *env,
                                                 jobject jobj,
                                                 jlong native_pointer,
                                                 jint j_length) {
    if (NULL == native_pointer) {
        LOGE("YOU MUST INIT BEFER USE");
        return NULL;
    }

    CDelayBuf *delay_buffer = reinterpret_cast<CDelayBuf *>(native_pointer);
    jbyteArray j_out_array = env->NewByteArray(j_length);
    jbyte *data = env->GetByteArrayElements(j_out_array,JNI_FALSE);

    int result = delay_buffer->DelayBufGet(reinterpret_cast<int16_t *>(data), j_length/2);
    if (result == -1) {
        LOGE("GET ERROR");
        return NULL;
    }
    env->SetByteArrayRegion(j_out_array,0,j_length,data);
    return j_out_array;
}

extern "C"
JNIEXPORT jint JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_size(JNIEnv *env,
                                                  jobject jobj,
                                                  jlong native_pointer) {
    if (NULL == native_pointer) {
        LOGE("YOU MUST INIT BEFER USE");
        return -1;
    }

    return ((CDelayBuf *) native_pointer)->DelayBufedSize();
}


extern "C"
JNIEXPORT void JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_reset(JNIEnv *env,
                                                   jobject jobj,
                                                   jlong native_pointer) {
    if (NULL == native_pointer) {
        LOGE("YOU MUST INIT BEFER USE");
        return;
    }

    ((CDelayBuf *) native_pointer)->DelayBufReset();
}


extern "C"
JNIEXPORT void JNICALL
Java_com_sabinetek_swiss_jni_buffer_JitterBuffer_delete(JNIEnv *env,
                                                    jobject jobj,
                                                    jlong native_pointer) {
    if (native_pointer != NULL) {
        delete (CDelayBuf*)native_pointer;
    }
}
