package com.codec;

import android.util.Log;

import static com.sabinetek.swiss.jni.INative.CodecType;


/**
 * <p>jitterbuffer功能的对外使用工具。</p>
 * <p>单例实现</p>
 */
public class JitterBuffer {
    private final static int JITTER_MAX_DELAY = 200;
    //jitterbuffer create 参数
    private final static int FRAME_TIME_SBC = 16;
    private final static int FRAME_TIME_CELT = 20;
    //get data len
    private final static int FRAME_LEN_SBC = 3072;
    private final static int FRAME_LEN_CELT = 3840;
    //get delay time
    private final static long GET_DELAY_TIME_SBC = 17415000;//ns 17.415ms  3072 = 48000/1000 * 16 * 2 * 2  3072 = 44100/1000 * 17.4150 * 2 * 2
    private final static long GET_DELAY_TIME_CELT = 20000000;//ns 20ms     3840 = 48000/1000 * 20 * 2 * 2  3840 = 44100/1000 * 21.7687 * 2 * 2
//    private final static long GET_DELAY_TIME_CELT = 21768700;//ns 20ms     3840 = 48000/1000 * 20 * 2 * 2  3840 = 44100/1000 * 21.7687 * 2 * 2

    private int frame_time = FRAME_TIME_SBC;
    private int frame_len = FRAME_LEN_SBC;
    public static long GET_DELAY_TIME = GET_DELAY_TIME_SBC;

    private com.sabinetek.swiss.jni.buffer.JitterBuffer jitterBuffer;
    private boolean openJitterBuffer = false;
    private boolean isPutFull = false;
    private int putTimes = 0;

    public static JitterBuffer getInstance() {
        return JitterBuffer.BufferHolder.INSTANCE;
    }

    /**
     * 根据coderType构建jitterbuffer
     *
     * @param coderType 解码类型
     */
    public void createBuffer(CodecType coderType) {
        if (openJitterBuffer) {
            if (null != jitterBuffer) {
                jitterBuffer.destory();
            }
            Log.e("JitterProvider", "createBuffer:" + coderType);
            jitterBuffer = new com.sabinetek.swiss.jni.buffer.JitterBuffer();
            if (coderType == CodecType.SBC) {
                frame_time = FRAME_TIME_SBC;
                frame_len = FRAME_LEN_SBC;
                GET_DELAY_TIME = GET_DELAY_TIME_SBC;
            } else if (coderType == CodecType.OPUS) {
                frame_time = FRAME_TIME_CELT;
                frame_len = FRAME_LEN_CELT;
                GET_DELAY_TIME = GET_DELAY_TIME_CELT;
            }
            jitterBuffer.createBuffer(48000, 2, frame_time, JITTER_MAX_DELAY);
            isPutFull = false;
            putTimes = 0;
        }
    }

    /**
     * 添加数据
     *
     * @param pcm 待加入数据
     */
    public void put(byte[] pcm) {
        if (null != jitterBuffer && null != pcm && pcm.length > 0) {
            int Pl = jitterBuffer.put(pcm, pcm.length);
            if (0 == Pl){
                putTimes ++;
            }
            if (!isPutFull && putTimes >= 6){// false --> true
                Log.e("JitterProvider","putTimes = "+putTimes);
                isPutFull = true;
            }
        }
    }

    /**
     * 获取数据
     *
     * @return 获取到数据
     */
    public byte[] get() {
        if (null != jitterBuffer && isPutFull) {
            return jitterBuffer.get(frame_len);
        }
        return null;
    }

    /**
     * 关闭jitterbuffer
     */
    public void release() {
        if (null != jitterBuffer) {
            jitterBuffer.destory();
            jitterBuffer = null;
            openJitterBuffer = false;
            isPutFull = false;
            putTimes = 0;
        }
    }

    /**
     * 设置是否开启jitterbuffer
     *
     * @param openJitterBuffer 开启标识
     */
    public void setOpenJitterBuffer(boolean openJitterBuffer) {
        this.openJitterBuffer = openJitterBuffer;
    }

    /**
     * jitterbuffer是否开启
     *
     * @return openJitterBuffer 是否开启
     */
    public boolean isOpen() {
        return openJitterBuffer;
    }

    private static class BufferHolder {
        protected static JitterBuffer INSTANCE = new JitterBuffer();
    }
}
