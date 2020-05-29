package com.sabinetek.swiss.jni.codec.raw;

import com.sabinetek.swiss.jni.codec.CodecSbc;

/**
 * sbc 编解码器功能的封装
 */
public class CoderRawSBC implements ICoderRaw {

    private CodecSbc codecSbc;
    private byte[] data;

    @Override
    public void initData() {
        codecSbc = new CodecSbc();
        codecSbc.SbcInit(0);
        data = new byte[1024 * 12];

    }

    @Override
    public int initEncode(int sampleRate, int channels, int bitrate,
                          int complexity, int bitsPerSample) {
        return 0;
    }

    @Override
    public int initDecode(int sampleRate, int channels, int version) {
        return 0;
    }


    @Override
    public byte[] decode(byte[] inputArray) {
        if (codecSbc != null && inputArray != null && data != null) {
            int len = codecSbc.SbcDecode(inputArray, inputArray.length, data);
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
        if (codecSbc != null && inputArray != null && data != null) {
            int len = codecSbc.SbcEncode(inputArray, inputArray.length, data);
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
        if (codecSbc != null) {
            codecSbc.SbcExit();
            codecSbc = null;
        }
    }

    @Override
    public void closeEncode() {
        if (codecSbc != null) {
            codecSbc.SbcExit();
            codecSbc = null;
        }
    }
}
