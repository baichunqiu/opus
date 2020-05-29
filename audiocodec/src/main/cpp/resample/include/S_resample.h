#include "Mixer_Def.h"
#include "Mixer_Resample.h"
#include "Mixer_RS_Linear.h"
#include "Mixer_RS_KW.h"
#include <stdint.h>
#include<iostream>

#ifndef S_RESAMPLE_H
#define S_RESAMPLE_H

typedef struct LIBPCM_RESAMPLE_EXE_SRC_PARAM_ST
{
	void * src_data;
	int32_t src_data_len;			//src data length
	int32_t src_used_data_len;	//src used data length
}LIBPCM_RESAMPLE_EXE_SRC_PARAM_T;

typedef struct LIBPCM_RESAMPLE_EXE_DES_PARAM_ST
{
	short * des_data;
	int32_t des_data_len;			//des data length
	int32_t des_used_data_len;	//des used data length
}LIBPCM_RESAMPLE_EXE_DES_PARAM_T;

class S_resample{
private:
	MIXER_RS_KW* m_kw_resample;
	MIXER_RS_LINEAR* m_liner_resample;
	uint32_t m_nBit;
	uint32_t m_nChannel;
	uint32_t m_nNewSampleRate;
	uint32_t m_nOldSampleRate;
	int32_t m_nLengthSrc;
public:
//	S_resample(uint32_t nNewSampleRate=44100,uint32_t nOldSampleRate=48000,uint32_t nBit = 16,uint32_t nChannel = 2,uint32_t nLengthSrc= 0);
	S_resample(uint32_t nNewSampleRate,uint32_t nOldSampleRate,uint32_t nBit,uint32_t nChannel,uint32_t nLengthSrc);
	~S_resample(){delete m_kw_resample;delete m_liner_resample;}
	bool s_resample_reset();
	bool s_resample_ioctl();
    bool s_resample_execute(LIBPCM_RESAMPLE_EXE_SRC_PARAM_T &src_param,  LIBPCM_RESAMPLE_EXE_DES_PARAM_T &des_param);
	
};

#endif
