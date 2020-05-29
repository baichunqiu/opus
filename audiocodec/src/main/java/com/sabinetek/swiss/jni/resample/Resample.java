package com.sabinetek.swiss.jni.resample;

/**
 * 重采样功能封装集
 */
public class Resample {

    static {
        System.loadLibrary("jni-resample");
    }

    private long nativePointer;

    /**
     * 构建重采样对象
     * @param srcSampleRate 原采样率
     * @param desSampleRate 目标采样率
     * @param channel 数据的声道
     */
    public Resample(int srcSampleRate,
                    int desSampleRate,
                    int channel) {
        nativePointer = init(srcSampleRate, desSampleRate, channel);
    }

    /**
     * 重采样
     * @param src 带采样数据
     * @param len 数据长度
     * @return 重采样数据
     */
    public byte[] resample(byte[] src,
                           int len) {
        return resample(nativePointer, src, len);
    }

    /**
     * 关闭重采样
     */
    public void close() {
        if (nativePointer != 0)
            close(nativePointer);
    }


    private native long init(int srcSampleRate,
                             int desSampleRate,
                             int channel);

    private native byte[] resample(long objPointer, byte[] src,
                                   int len);

    private native void close(long objPointer);

}
