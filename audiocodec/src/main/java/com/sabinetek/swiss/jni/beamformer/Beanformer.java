package com.sabinetek.swiss.jni.beamformer;

/**
 * <p>波束形成功能封装。</p>
 * <p>目前主要使用转单声道功能</p>
 */
public class Beanformer {

    private boolean isClose = true;

    // 1 frame size 960 samples
    private static final int FRAME_SIZE = 960;

    // first 2 is short convert byte, second 2 is channel
    private static final int FRAME_SIZE_IN_BYTES = FRAME_SIZE * 2 * 2;

    private long nativePointer = 0;

    static {
        System.loadLibrary("jni-beamformer");
    }

    /**
     * 构建波束形成器实例
     */
    public Beanformer() {
        nativePointer = open(44100, FRAME_SIZE);
        isClose = false;
    }

    /**
     * 波束形成器 数据处理
     *
     * @param pcm 需要处理的PCM数据
     * @param len 需要处理的PCM数据长度
     * @return 处理后的数
     */
    public byte[] process(byte[] pcm, int len) {
        if (isClose)
            return null;

        int frames = len / FRAME_SIZE_IN_BYTES;

        // if one frame
        if (frames == 1) {
            return nativeProcess(nativePointer, pcm, FRAME_SIZE_IN_BYTES);
        } else {
            //  Divide by 2 is convert mono
            byte[] totalFrames = new byte[FRAME_SIZE_IN_BYTES * frames / 2];
            for (int i = 0; i < frames; i++) {
                byte[] bytes = new byte[FRAME_SIZE_IN_BYTES];
                System.arraycopy(pcm,
                        i * FRAME_SIZE_IN_BYTES,
                        bytes,
                        0,
                        FRAME_SIZE_IN_BYTES);

                byte[] processedPcm = nativeProcess(nativePointer, bytes, FRAME_SIZE_IN_BYTES);

                System.arraycopy(processedPcm,
                        0,
                        totalFrames,
                        i * processedPcm.length,
                        processedPcm.length);
            }

            return totalFrames;
        }
    }

    /**
     * 是否波束形成器
     */
    public void release() {
        if (nativePointer != 0) {
            nativeClose(nativePointer);
            isClose = true;
        }
    }

    private native long open(int sampleRate, int frameSize);

    private native byte[] nativeProcess(long nativePointer, byte[] inPcm, int len);

    private native void nativeClose(long nativePointer);
}
