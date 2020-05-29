
#include <jni.h>
#include <stdio.h>

#include "Def_DataType.h"
#include "S_resample.h"
#include "../../../../../../libmedia/src/main/cpp/base_include.h"
#include "../include/S_resample.h"

#define TAG    "resample"


S_resample *cS_resample = NULL;
int _outSR = 0, _inSR = 0;

extern "C"
JNIEXPORT jlong JNICALL
Java_com_sabinetek_swiss_jni_resample_Resample_init(JNIEnv *env,
                                                jobject jobj,
                                                jint inSR,
                                                jint outSR,
                                                jint channel) {

    S_resample *cS_resample = new S_resample(outSR, inSR, 16, channel, 0);
    return reinterpret_cast<jlong>(cS_resample);
}

extern "C"
JNIEXPORT jbyteArray JNICALL
Java_com_sabinetek_swiss_jni_resample_Resample_resample(JNIEnv *env,
                                                    jobject jobj,
                                                    jlong nativePointer,
                                                    jbyteArray inparam,
                                                    jint inLength) {


    jbyte *ins = env->GetByteArrayElements(inparam, JNI_FALSE);
    jbyteArray out = env->NewByteArray(1024 * 5);
    jbyte *outs = env->GetByteArrayElements(out, JNI_FALSE);


    LIBPCM_RESAMPLE_EXE_SRC_PARAM_ST strInbuffe;
    strInbuffe.src_data_len = inLength; //src_data_len;
    strInbuffe.src_used_data_len = 0;
    strInbuffe.src_data = (jshort *) ins;

    LIBPCM_RESAMPLE_EXE_DES_PARAM_ST strOutbuffe;
    strOutbuffe.des_data_len = 1024 * 5; //des_data_len;
    strOutbuffe.des_used_data_len = 0;
    strOutbuffe.des_data = (jshort *) outs;

    reinterpret_cast<S_resample *>(nativePointer)->s_resample_execute(strInbuffe, strOutbuffe);

    jbyteArray result = env->NewByteArray(strOutbuffe.des_used_data_len);

    env->SetByteArrayRegion(result,
                            0,
                            strOutbuffe.des_used_data_len,
                            (const jbyte *) strOutbuffe.des_data);

    env->ReleaseByteArrayElements(inparam, ins, 0);
    env->ReleaseByteArrayElements(out, outs, 0);
    env->DeleteLocalRef(out);
    return result;
}


extern "C"
JNIEXPORT void JNICALL
Java_com_sabinetek_swiss_jni_resample_Resample_close__J(JNIEnv *env,
                                                 jobject jobj,
                                                 jlong nativePointer) {
    if (nativePointer != 0) {
        S_resample *s_resample = reinterpret_cast<S_resample *>(nativePointer);
        delete s_resample;
    }
}




