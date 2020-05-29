package com.sabinetek.swiss.jni.codec;

/**
 * celt 编解码器
 */
public class CodecOpus {


    static {
        System.loadLibrary("jni-opus");
    }

    public static native int initEncode(int sampleRate, int channels, int bitrate,
                                        int complexity, int bitsPerSample);

    public static native int initDecode(int sampleRate, int channels,int version);

    public static native int encode(byte[] inputArray, int frameSize,
                                    byte[] outputArray);

    public static native int decode(byte[] inputArray, int frameSize,
                                    byte[] outputArray);

    public static native void closeEncode();

    public static native void closeDecode();


}
