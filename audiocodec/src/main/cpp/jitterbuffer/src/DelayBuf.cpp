#include "DelayBuf.h"


CDelayBuf::CDelayBuf() :m_nSamplingRate(48000), m_nChannel(2), m_nFrameTime(20), m_nMaxDelayTime(200)
{
	m_nMaxLevel = 0;
	m_nInstLevel = 0;
	m_lastOp = OP_GET;

	m_nRecalcLevelTimer = 0;
	m_nEffCnt = 0;

	init();
}


CDelayBuf::~CDelayBuf()
{
	delete m_pWsola;
	delete m_pCircBuf;
	delete m_pShrinkBuf;
	delete m_pExpandBuf;
#ifdef SHOWLOG
	delete m_locallog;
#endif
}
CDelayBuf::CDelayBuf(uint32_t nSamplingRate,
	uint32_t	nChannel,
	uint32_t	nFrameTime,//in ms
	uint32_t	nMaxDelay//in ms
	)
{
	m_nSamplingRate = nSamplingRate;
	m_nChannel = nChannel;
	m_nFrameTime = nFrameTime;
	m_nMaxDelayTime = nMaxDelay;

	m_nFadeEffectTime = 16;

	m_nMaxLevel = 0;
	m_nInstLevel = 0;
	m_lastOp = OP_GET;

	m_nRecalcLevelTimer = 0;
	m_nEffCnt = 0;

	init();
}
void CDelayBuf::init()
{
	m_nSamplesPerFrame = m_nSamplingRate / 1000 * m_nFrameTime * m_nChannel;
	m_nMaxDelaySample = m_nSamplingRate / 1000 * m_nMaxDelayTime * m_nChannel;

	m_nExpandBufSize = m_nSamplesPerFrame * 6;

	m_nShrinkBufSize = m_nMaxDelaySample * 4;
	m_nCircBufSize = m_nMaxDelaySample * 4;

	m_nMinimumExpandLen = m_nSamplesPerFrame * 2;

	m_nFadeEffectSample = m_nSamplingRate / 1000 * m_nFadeEffectTime * m_nChannel;

	m_pWsola = new CWsola(m_nSamplingRate, m_nChannel, m_nFrameTime, 5);
	m_pCircBuf = new CCircBuffer(m_nCircBufSize);

	m_pShrinkBuf = new int16_t[m_nShrinkBufSize];
	m_pExpandBuf = new int16_t[m_nExpandBufSize];

	m_mutex =  (pthread_mutex_t)PTHREAD_MUTEX_INITIALIZER;

#ifdef SHOWLOG
	m_locallog = new CLocalLog("DelayBuf.log");
#endif
}
int16_t CDelayBuf::fade_out(int16_t* pcm_buf, uint32_t nLen)
{
	if (nLen < m_nFadeEffectSample)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf:  fade out can not reach here");
#endif
		return -1;
	}
	//do fade out
	uint32_t	nIndex = 0;
	int16_t		nData = 0;
	for (nIndex = 0; nIndex < m_nFadeEffectSample; nIndex++)
	{
		nData = *(pcm_buf + nLen - m_nFadeEffectSample + nIndex);
		nData = nData * ((float)(m_nFadeEffectSample - nIndex) / (float)m_nFadeEffectSample);
		*(pcm_buf + nLen - m_nFadeEffectSample + nIndex) = nData;
	}

	return 0;
}
int16_t CDelayBuf::fade_in(int16_t* pcm_buf, uint32_t nLen)
{
	if (nLen < m_nFadeEffectSample)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf:  fade in can not reach here");
#endif
		return -1;
	}
	//do fade out
	uint32_t	nIndex = 0;
	int16_t		nData = 0;
	for (nIndex = 0; nIndex < m_nFadeEffectSample; nIndex++)
	{
		nData = *(pcm_buf + nIndex);
		nData = nData * ((float)(nIndex) / (float)m_nFadeEffectSample);
		*(pcm_buf + nIndex) = nData;
	}

	return 0;
}

void CDelayBuf::update(enum OP op)
{
	uint32_t nEraseCnt = 0;
	uint32_t nErasedCnt = 0;
	uint32_t nExpandCnt = 0;
	uint32_t nExpandedCnt = 0;
	uint32_t nPCMLen = 0;
	int32_t	nNewEffCnt = 0;

	
	m_nRecalcLevelTimer -= (m_nInstLevel * m_nFrameTime) >> 1;

	if (op == m_lastOp)
	{
		m_nInstLevel++;
	}
	else
	{
		m_nInstLevel = 1;
	}
	m_lastOp = op;
	//According to the level,shrink the buffer,or expand buffer
	if (m_nInstLevel > m_nMaxLevel)
	{
		m_nMaxLevel = m_nInstLevel;
	}

	if (m_nRecalcLevelTimer < 0)
	{
		//We should recalculate the timer counter
		//At least we should have 2 frames for shrink or expand
		nNewEffCnt = (m_nMaxLevel + 2)* m_nSamplesPerFrame;
		SMOOTH(m_nEffCnt, nNewEffCnt); //Smooth the m_nEffCnt
		if (m_nEffCnt % m_nChannel)
			m_nEffCnt = nNewEffCnt + 1;
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf:  recalc effcnt %d level %d/%d", nNewEffCnt, m_nInstLevel, m_nMaxLevel);
#endif
		m_nMaxLevel = 0;
		m_nRecalcLevelTimer = 2000;//We do the calc every 2000ms
	}

	if (op == OP_GET)
	{
		nPCMLen = m_pCircBuf->CircBufGetLen();
		if (nPCMLen <= m_nSamplesPerFrame)
		{
#ifdef SHOWLOG
			m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: data is too small,we cannot perform expand");
#endif
		}
		else if ((nPCMLen - m_nSamplesPerFrame) <= m_nMinimumExpandLen)
		{
#ifdef SHOWLOG
			m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: data is critical,so we need expand");
#endif
			nExpandCnt = m_nSamplesPerFrame;
			nExpandedCnt = expand_buffer(nExpandCnt);
#ifdef SHOWLOG
			m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: expand %d/%d", nExpandedCnt, nExpandCnt);
#endif
			//check again
			nPCMLen = m_pCircBuf->CircBufGetLen();
			if (nPCMLen <= m_nMinimumExpandLen)
			{
#ifdef SHOWLOG
				m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: expand fail,this will cause pop");
#endif
				nExpandCnt = m_nSamplesPerFrame;
				m_pCircBuf->CircBufAdvWritePtr(nExpandCnt);
			}
			m_nRecalcLevelTimer = -1;
			m_lastGetFrameOP = OP_EXPAND_DATA;
		}
		else
		{
			if (m_lastGetFrameOP == OP_EXPAND_DATA)
			{
				//we should do a fade in the put in data
			}
			m_lastGetFrameOP = OP_NORMAL_DATA;
		}
	}
	else if (op == OP_PUT)
	{
		nPCMLen = m_pCircBuf->CircBufGetLen();

		if (nPCMLen > (m_nSamplesPerFrame + m_nEffCnt))
		{
			nEraseCnt = m_nSamplesPerFrame >> 1;

			if (nEraseCnt % m_nChannel != 0)
				nEraseCnt += 1;

			if (nPCMLen < m_nMaxDelaySample)
			{
#ifdef SHOWLOG
				m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: no need shrink");
#endif
			}
			else
			{
				nErasedCnt = shrink_buffer(nEraseCnt);
#ifdef SHOWLOG
				m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: shrink %d/%d buflen %d effect %d", nErasedCnt, nEraseCnt, nPCMLen, m_nEffCnt);
#endif
			}

			//check again
/*			nPCMLen = m_pCircBuf->CircBufGetLen();
			if (nPCMLen + m_nSamplesPerFrame > m_nMaxDelaySample)
			{
				//earse is fail
				nEraseCnt = nPCMLen + m_nSamplesPerFrame - m_nMaxDelaySample;
				m_pCircBuf->CircBufAdvReadPtr(nEraseCnt);
#ifdef SHOWLOG
				m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: erase fail,adv read ptr");
#endif
			}*/
		}
	}
}
int16_t CDelayBuf::shrink_buffer(uint32_t erase_cnt)
{
	int16_t *buf1, *buf2;
	uint32_t buf1len = 0;
	uint32_t buf2len = 0;
	uint32_t nPCMLen = 0;
	uint32_t nErasedCnt = 0;

	if (m_pCircBuf == NULL || m_pWsola == NULL)
		return 0;
	
	nPCMLen = m_pCircBuf->CircBufGetLen();
	if (nPCMLen == 0)
		return 0;

	m_pCircBuf->CircBufPeek(m_pShrinkBuf, nPCMLen);
	nErasedCnt = m_pWsola->wsola_discard(m_pShrinkBuf, nPCMLen, m_nShrinkBufSize, erase_cnt);
	m_pCircBuf->CircBufSetLen(0);
	m_pCircBuf->CircBufWrite(m_pShrinkBuf, nPCMLen - nErasedCnt);

	return nErasedCnt;
}


int16_t CDelayBuf::expand_buffer(uint32_t expand_cnt)
{
	int16_t *buf1, *buf2;
	uint32_t buf1len = 0;
	uint32_t buf2len = 0;
	uint32_t nPCMLen = 0;
	uint32_t nExpandedCnt = 0;

	if (m_pCircBuf == NULL || m_pWsola == NULL)
		return 0;

	nPCMLen = m_pCircBuf->CircBufGetLen();
	if (nPCMLen == 0)
		return 0;

	m_pCircBuf->CircBufPeek(m_pExpandBuf, nPCMLen);
	nExpandedCnt = m_pWsola->wsola_generate(m_pExpandBuf, nPCMLen, m_nExpandBufSize, m_nSamplesPerFrame);
	if (m_lastGetFrameOP == OP_EXPAND_DATA)
	{
		//do a fade out to generated data
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: Expand Fadeout %d", nExpandedCnt);
#endif
		fade_out(m_pExpandBuf, nPCMLen + nExpandedCnt);
	}
	else
	{
		//do nothing
	}

	m_pCircBuf->CircBufSetLen(0);
	m_pCircBuf->CircBufWrite(m_pExpandBuf, nPCMLen + nExpandedCnt);

	return nExpandedCnt;
}

int16_t CDelayBuf::DelayBufPut(int16_t* pcm_buf, uint32_t nLen)
{
#ifdef SHOWLOG
	m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: Put");
#endif
	if (nLen != m_nSamplesPerFrame)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: samples is not equal");
#endif
		return -1;
	}
	pthread_mutex_lock(&m_mutex);
	update(OP_PUT);
	uint32_t nPCMLen = 0;
	nPCMLen = m_pCircBuf->CircBufGetLen();
	/*Overflow Checking*/
	if (nPCMLen + m_nSamplesPerFrame > m_nCircBufSize)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: no enough space,this can not be true");
#endif
	}
	else
	{
		if (m_lastGetFrameOP == OP_EXPAND_DATA)
		{
			//this should do a face in data
			fade_in(pcm_buf, nLen);
			m_lastGetFrameOP = OP_NORMAL_DATA;
		}
		m_pCircBuf->CircBufWrite(pcm_buf, nLen);
	}
	pthread_mutex_unlock(&m_mutex);

	return 0;
}
int16_t CDelayBuf::DelayBufGet(int16_t* pcm_buf, uint32_t nLen)
{
#ifdef SHOWLOG
	m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: Get");
#endif
	if (nLen != m_nSamplesPerFrame)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: get samples is not equal");
#endif
		return m_nSamplesPerFrame;
	}
	uint32_t nPCMLen = 0;

	pthread_mutex_lock(&m_mutex);

	update(OP_GET);//When perform op_get,we may perform the expand operation

	nPCMLen = m_pCircBuf->CircBufGetLen();
	if (nPCMLen <= m_nSamplesPerFrame)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Delaybuf: no enough data,this can not be true");
#endif
		pthread_mutex_unlock(&m_mutex);
		return -1;
	}
	else
		m_pCircBuf->CircBufRead(pcm_buf, nLen);

	pthread_mutex_unlock(&m_mutex);

	return 0;
}

uint32_t CDelayBuf::DelayBufedSize()
{
	uint32_t nPCMLen = 0;

	pthread_mutex_lock(&m_mutex);
	nPCMLen = m_pCircBuf->CircBufGetLen();
	pthread_mutex_unlock(&m_mutex);

	return nPCMLen / m_nSamplesPerFrame;
}
int16_t CDelayBuf::DelayBufReset()
{
	m_pCircBuf->CircBufReset();
	return 0;
}
