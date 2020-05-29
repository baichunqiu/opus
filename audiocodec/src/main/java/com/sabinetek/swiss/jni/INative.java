package com.sabinetek.swiss.jni;

/**
 * Native 层接口
 */
public interface INative {
    /**
     * 构建codec
     * @param coderType 枚举类型CodecType
     * @param sampleRate 采样率
     * @param channels 声道数
     * @param version 版本号
     * @return 是否构建新的codec
     */
    boolean buildCodec(CodecType coderType, int sampleRate, int channels, int version);

    /**
     * 解码
     * @param dArr 待解码数据
     * @return 解码后数据
     */
     byte[] decode(byte[] dArr);

    /**
     * 编码
     * @param dArr
     * @return
     */
     byte[] encode(byte[] dArr);

    /**
     * 释放解码器
     */
     void releaseCodec();

    /**
     * 音频PCM数据转换 实现机制：采样率非默认则重采样，声道数非默认则使用波束形成器转声道
     *
     * @param orgPcm  待处理的原始数据
     * @param sample  目标采样率
     * @param channel 目标声道数
     * @return 转换后PCM
     */
    byte[] conversion(byte[] orgPcm, int sample, int channel);

    /**
     * 释放转换器
     */
    void releaseConversion();

    /**
     * 波束形成
     *
     * @param pcm org data 待波形形成的原始数据
     * @return 波束形成数据
     */
    byte[] beanformer(byte[] pcm);

    /**
     * 释放波束形成器
     */
    void releaseBeamformer();

    /**
     * 重采样（只修改采样率，声道数不变）
     *
     * @param pcm       待重采样原始数据
     * @param oldsample 原采样率
     * @param newsample 目标采样率
     * @param channle   数据的声道
     * @return 重采样数组
     */
    byte[] resample(byte[] pcm, int oldsample, int newsample, int channle);

    /**
     * 释放重采样
     */
    void releaseResample();

    /**
     * 解码器类型枚举
     */
    enum CodecType {
        SBC,
        OPUS
    }
}
