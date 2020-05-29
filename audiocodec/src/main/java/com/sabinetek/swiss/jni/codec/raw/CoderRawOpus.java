package com.sabinetek.swiss.jni.codec.raw;

import com.sabinetek.swiss.jni.codec.CodecOpus;

/**
 * opus 编解码器功能的封装
 */
public class CoderRawOpus implements ICoderRaw {
    private CodecOpus codecOpus;
    private byte[] data;

    @Override
    public void initData() {
        codecOpus = new CodecOpus();
        data = new byte[1024 * 12];
    }

    @Override
    public int initEncode(int sampleRate, int channels, int bitrate,
                          int complexity, int bitsPerSample) {
        //0 表示初始化成功 -1 表示初始化失败
        int state = -1;
        if (codecOpus != null) {
            state = codecOpus.initEncode(sampleRate, channels, bitrate, complexity, bitsPerSample);
        }
        return state;
    }

    @Override
    public int initDecode(int sampleRate, int channels, int version) {
        //0 表示初始化成功 -1 表示初始化失败
        int state = -1;
        if (codecOpus != null) {
            state = codecOpus.initDecode(sampleRate, channels, version);
        }
        return state;
    }

    @Override
    public byte[] decode(byte[] inputArray) {
        if (codecOpus != null && inputArray != null && data != null) {
            int len = codecOpus.decode(inputArray, inputArray.length, data);
            if (len > 0) {
                inputArray = new byte[len];
                System.arraycopy(data, 0, inputArray, 0, inputArray.length);
                return inputArray;
            }
        }
        return null;
    }

    @Override
    public byte[] encode(byte[] inputArray) {
        if (codecOpus != null && inputArray != null && data != null) {
            int len = codecOpus.encode(inputArray, inputArray.length, data);
            if (len > 0) {
                inputArray = new byte[len];
                System.arraycopy(data, 0, inputArray, 0, inputArray.length);
                return inputArray;
            }
        }
        return null;
    }

    @Override
    public void closeDecode() {
        if (codecOpus != null) {
            codecOpus.closeDecode();
        }
    }

    @Override
    public void closeEncode() {
        if (codecOpus != null) {
            codecOpus.closeEncode();
        }
    }
}
