#ifndef _DELAY_BUF_H
#define _DELAY_BUF_H

#include "datatypes.h"
#ifdef SHOWLOG
#include "locallog.h"
#endif
#include "Wsola.h"
#include "CircBuffer.h"

#define SMOOTH_UP(cur, target) cur = (cur + target*3) >> 2
#define SMOOTH_DOWN(cur, target) cur = (cur*3 + target) >> 2
#define SMOOTH(cur, target) \
    if (cur < target) SMOOTH_UP(cur, target); \
	else SMOOTH_DOWN(cur, target)


class CDelayBuf
{
public:
	enum OP
	{
		OP_PUT,
		OP_GET
	};
	enum WSOLA_OP
	{
		OP_EXPAND_DATA,
		OP_SHRINK_DATA,
		OP_NORMAL_DATA
	};
public:
	CDelayBuf(uint32_t nSamplingRate,
		uint32_t	nChannel,
		uint32_t	nFrameTime,//in ms
		uint32_t	nMaxDelay//in ms
		);
	CDelayBuf();
	~CDelayBuf();
private:
	void init();
	void update(enum OP op);
	int16_t shrink_buffer(uint32_t erase_cnt);
	int16_t expand_buffer(uint32_t expand_cnt);


	int16_t fade_out(int16_t* pcm_buf,uint32_t nLen);
	int16_t fade_in(int16_t* pcm_buf, uint32_t nLen);

public:
	int16_t DelayBufPut(int16_t* pcm_buf, uint32_t nLen);
	int16_t DelayBufGet(int16_t* pcm_buf, uint32_t nLen);
	uint32_t DelayBufedSize();
	int16_t DelayBufReset();

private:
	CWsola*			m_pWsola;
	CCircBuffer*	m_pCircBuf;
	int16_t*		m_pShrinkBuf;
	int16_t*		m_pExpandBuf;

private:
	uint32_t	m_nSamplingRate;
	uint32_t	m_nChannel;
	uint32_t	m_nFrameTime;
	uint32_t	m_nMaxDelayTime;
	uint32_t	m_nFadeEffectTime;

private:
	uint32_t	m_nSamplesPerFrame;
	uint32_t	m_nMaxDelaySample;

	uint32_t	m_nMinimumExpandLen;

	uint32_t	m_nExpandBufSize;
	uint32_t	m_nShrinkBufSize;

	uint32_t	m_nCircBufSize;

	uint32_t	m_nFadeEffectSample;

private:

	uint32_t	m_nMaxLevel;	//Maxlevel of the data,this is smoothed
	uint32_t	m_nInstLevel;	//Instantaneous level of current data
	OP			m_lastOp;		//Last operation of put or get
	WSOLA_OP	m_lastGetFrameOP;	//Last WSola op for fade in /out operation
	WSOLA_OP	m_lastPutFrameOP;	//Last WSola op for fade in /out operation

	int32_t		m_nRecalcLevelTimer;
	int32_t		m_nEffCnt;

private:
    mutable pthread_mutex_t m_mutex;

private:
#ifdef SHOWLOG
	CLocalLog	*m_locallog;
#endif
};

#endif
