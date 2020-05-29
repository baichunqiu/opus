package com.sabinetek.swiss.jni.buffer;

/**
 * <p>JitterBuffer 功能的封装</p>
 * <p>功能：实现减小数据流的波动，接近恒定的时间间隔获取swiss的pcm数据流。</p>
 */
public class JitterBuffer {

    static {
        System.loadLibrary("jni-jitterbuffer");
    }

    private long mNativePointer;

    /**
     * 构建jitterBuffer
     *
     * @param srcSampleRate 数据采样率
     * @param channel       数据声道数
     * @param frameTime     帧时长
     * @param nMaxDelay     最大延迟时间
     */
    public void createBuffer(int srcSampleRate,
                             int channel,
                             int frameTime,
                             int nMaxDelay) {
        mNativePointer = init(srcSampleRate, channel, frameTime, nMaxDelay);
    }

    /**
     * 向jitterBuffer 添加数据
     *
     * @param data   添加PCM数据 byte[]
     * @param length 长度
     * @return 0 标识成功；-1 失败。
     */
    public int put(byte[] data,
                   int length) {
        return put(mNativePointer, data, length);
    }

    /**
     * 从jitterBuffer取出数据
     *
     * @param length 帧时长
     * @return byte[] 可能为null
     */
    public byte[] get(int length) {
        return get(mNativePointer, length);
    }

    /**
     * 获取jitterbuffer的缓存大小
     *
     * @return jitterbuffer当前的缓存数据大小
     */
    public int size() {
        return size(mNativePointer);
    }

    /**
     * 重置buffer
     */
    public void reset() {
        reset(mNativePointer);
    }

    /**
     * 销毁jitterbuffer
     */
    public void destory() {
        delete(mNativePointer);
    }

    private native long init(int srcSampleRate,
                             int channel,
                             int frameTime,
                             int nMaxDelay);

    private native int put(long nativePointer,
                           byte[] data,
                           int length);

    private native byte[] get(long nativePointer,
                              int length);

    private native int size(long nativePointer);

    private native void reset(long nativePointer);

    private native void delete(long nativePointer);
}
