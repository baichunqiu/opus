package com.sabinetek.swiss.jni.codec;

/**
 * sbc 编解码器
 */
public class CodecSbc {

    static {
        System.loadLibrary("jni-sbc");
    }

    public static native void SbcInit(int flag);

    public static native int SbcDecode(byte[] input, int inputlen,
                                       byte[] output);

    public static native int SbcEncode(byte[] input, int inputlen,
                                       byte[] output);

    public static native void SbcExit();
}
