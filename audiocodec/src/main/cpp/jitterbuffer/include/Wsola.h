#ifndef _WSOLA_H
#define _WSOLA_H

#include "datatypes.h"
#ifdef SHOWLOG
#include "locallog.h"
#endif

class CWsola
{
public:
	CWsola(uint16_t nSamplingRate,uint16_t nChannel,uint16_t nFramesInterval,uint16_t nWindow_Size);
	CWsola();
	~CWsola();
public:
	/*
	Description
		this function will generate some samples
	Input
		buf: the original pcm data,this buf should at least contain 1.5 frames of data
		nPCMLen: the buf pcm data length
		nBufLen: the buf total length,should at least 4 frames length
		nExpandedLen: the expand at least
	Return
		the expanded data length
	*/
	int16_t wsola_generate(int16_t buf[], uint32_t nPCMLen, uint32_t nBufLen, uint32_t nExpandedLen);
	/*
	Description
		this function will discard some samples
	Input
		buf: the original pcm data,this buf should at least contain 1.5 frames of data
		nPCMLen: the buf pcm data length
		nBufLen: the buf total length,should at least 4 frames length
		nCompressedLen: the discard at least
	Return
		the compressed data length
	*/
	int16_t wsola_discard(int16_t buf[], uint32_t nPCMLen, uint32_t nBufLen, uint32_t nCompressedLen);
private:
	/*
	Description
		this function will compress the buf for del_cnt
	Input
		pcm_buf: the original pcm data
		buf_size: the buf_size
		del_cnt: the compress size,in mono or stereo samples
	Return
		the deleted data length
	*/
	uint16_t compress(int16_t *pcm_buf, uint16_t pcm_size, uint16_t del_cnt);
	uint16_t compress(int16_t *pcm_buf_l, int16_t *pcm_buf_r, uint16_t pcm_size, uint16_t del_cnt);
	/*
	Description
		this function will expand the buf for needed data length
	Input
		pcm_buf: the original pcm data
		pcm_size: the buf_size
		needed: the expand size,in mono or stereo samples
	Return
		the expand data length
	*/
	uint16_t expand(int16_t *pcm_buf, uint16_t pcm_size, uint16_t needed);
	uint16_t expand(int16_t *pcm_buf_l, int16_t *pcm_buf_r,uint16_t pcm_size, uint16_t needed);

private:
	/*********************************************
	Description
		this function find the audio pitch
	Inputs
		from: the original pcm data, which is mono
		begin: from where to find the pitch
		end: to where will stop find the pitch
		template_size:the signal template size
		first: we find the local optimize or global optimized data
	Return
		The best start pos of the template data
	**********************************************/
	int16_t *find_pitch(int16_t *from, int16_t *begin, int16_t *end,
		uint16_t template_size, int16_t first);
	int16_t *find_pitch(int16_t *from, int16_t *begin, int16_t *end,
		uint16_t template_size, int16_t first, double* corr_value);
	/*********************************************
	Description
		this function will be used in expand and compress
	Inputs
		windows_size: the window size,which is half of the whold window
		templ:the original pcm data, which is mono
		signal_match:the finded max corr pos
		window: hanning window or other windows,while is only half of matlab window
	Outputs
		dst:the overlap add data
	Return
		null
	**********************************************/
	void overlapp_add(int16_t dst[], uint16_t window_size,
		int16_t templ[], int16_t signal_match[], float window[]);
	/*********************************************
	Description
		this function will create the hanning window
	Inputs
		pw:the hanning window data to be put in
		window_size: the hanning window size
	Outputs
		null
	Return
		null
	**********************************************/
	void create_win(float *pw, uint16_t window_size);
private:
	void init(); //allocate the buffer and init some buffer


private:
	float*	m_pHanning_Window;
	int16_t*	m_pMergeBuf;//This is used for window fade in/out
	int16_t*	m_pTempExpandL;	//This is a temp buffer
	int16_t*	m_pTempExpandR;	//This is a temp buffer
	int16_t*	m_pTempCompressL;	//This is a temp buffer
	int16_t*	m_pTempCompressR;	//This is a temp buffer

private:
	uint16_t m_nSamplingRate; 
	uint16_t m_nChannel;
	uint16_t m_nFramesInterval;
	uint16_t m_nWindow_Size;
private:
	//TemplateSize should always <= Hanning Window Size
	uint16_t m_nSamplesPerFrame;
	uint16_t m_nTemplateSize;
	uint16_t m_nExpand_Max_Dist;
	uint16_t m_nExpand_Min_Dist;
	uint16_t m_nTempBufSize;
private:
#ifdef SHOWLOG
	CLocalLog	*m_locallog;
#endif
};

#endif
