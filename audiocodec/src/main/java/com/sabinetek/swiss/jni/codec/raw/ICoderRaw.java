package com.sabinetek.swiss.jni.codec.raw;

/**
 * 编解码器接口
 */
public interface ICoderRaw {

    /**
     * 初始化数据
     */
    public void initData();

    /**
     * 初始化编码器
     * @param sampleRate 采样率
     * @param channels 声道数
     * @param bitrate 比特率
     * @param complexity 复杂性
     * @param bitsPerSample 字节预留
     * @return 0：成功 -1：失败
     */
    public int initEncode(int sampleRate, int channels, int bitrate,
                          int complexity, int bitsPerSample);

    /**
     * 初始化解码器
     * @param sampleRate 采样率
     * @param channels 声道数
     * @param version codec版本
     * @return 0：成功 -1：失败
     */
    public int initDecode(int sampleRate, int channels, int version);

    /**
     * 解码
     * @param inputArray 待解码数据
     * @return 解码后PCM数据
     */
    public byte[] decode(byte[] inputArray);

    /**
     * 编码
     * @param inputArray 待编码数据
     * @return 编码后数组
     */
    public byte[] encode(byte[] inputArray);

    /**
     * 关闭解码器
     */
    public void closeDecode();

    /**
     * 关闭编码器
     */
    public void closeEncode();

}
