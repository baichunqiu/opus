package com.codec;

import android.util.Log;

import com.sabinetek.swiss.jni.INative;
import com.sabinetek.swiss.jni.beamformer.Beanformer;
import com.sabinetek.swiss.jni.codec.raw.CoderRawOpus;
import com.sabinetek.swiss.jni.codec.raw.CoderRawSBC;
import com.sabinetek.swiss.jni.codec.raw.ICoderRaw;
import com.sabinetek.swiss.jni.resample.Resample;

/**
 * codec Provider
 */
public class CodecProvider implements INative {
    private final static String Tag = "NativeProvider";
    private ICoderRaw icoderRaw;
    private CodecType lastCodeType;
    private Beanformer beamformer;
    private Resample resample;
    public static int DEFAULT_SAMPLE = 44100;//默认采样率:44.1k
    public static int DEFAULT_CHANNEL = 2;//默认采声道：2

    @Override
    public boolean buildCodec(CodecType coderType, int sampleRate, int channels,  int bitrate) {
        if (null != icoderRaw) {
            if (coderType == lastCodeType) {
                return false;
            }
            Log.i(Tag, "MODIFY CODEC");
            releaseCodec();
        }
        switch (coderType) {
            case OPUS:
                icoderRaw = new CoderRawOpus();
                break;
            case SBC:
                icoderRaw = new CoderRawSBC();
                break;
            default:
        }
        if (null != icoderRaw) {
            icoderRaw.initData();
            int encode = icoderRaw.initEncode(sampleRate, channels, bitrate,
            1, 16);
            Log.i(Tag, "buildCodec:encode = " + encode);
            icoderRaw.initDecode(sampleRate, channels, 1);
        }
        lastCodeType = coderType;

        Log.i(Tag, "buildCodec:coderType = " + coderType);
        return true;
    }

    @Override
    public byte[] encode(byte[] dArr) {
        if (null != icoderRaw) {
            return icoderRaw.encode(dArr);
        }
        return null;
    }

    @Override
    public byte[] decode(byte[] dArr) {
        if (null != icoderRaw) {
            return icoderRaw.decode(dArr);
        }
        return null;
    }

    @Override
    public void releaseCodec() {
        Log.i(Tag, "releaseCodec");
        if (null != icoderRaw) {
            icoderRaw.closeDecode();
            icoderRaw.closeEncode();
            icoderRaw = null;
        }
    }

    @Override
    public byte[] conversion(byte[] orgPcm, int sample, int channel) {
        if (null == orgPcm) return null;
        if (DEFAULT_CHANNEL != channel) {
            orgPcm = beanformer(orgPcm);
        }
        if (DEFAULT_SAMPLE != sample) {
            orgPcm = resample(orgPcm, DEFAULT_SAMPLE, sample, channel);
        }
        return orgPcm;
    }

    @Override
    public void releaseConversion() {
        releaseBeamformer();
        releaseResample();
    }

    @Override
    public byte[] beanformer(byte[] pcm) {
        if (null == beamformer) {
            beamformer = new Beanformer();
        }
        if (null == pcm) return null;
        return beamformer.process(pcm, pcm.length);
    }

    @Override
    public void releaseBeamformer() {
        if (null != beamformer) {
            beamformer.release();
            beamformer = null;
        }
    }

    @Override
    public byte[] resample(byte[] pcm, int oldsample, int newsample, int channle) {
        if (null == resample) {
            resample = new Resample(oldsample, newsample, channle);
        }
        if (null == pcm) return null;
        return resample.resample(pcm, pcm.length);
    }

    @Override
    public void releaseResample() {
        if (null != resample) {
            resample.close();
            resample = null;
        }
    }
}
