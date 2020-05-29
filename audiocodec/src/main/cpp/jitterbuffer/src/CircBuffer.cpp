#include "CircBuffer.h"
#include <string.h>

CCircBuffer::CCircBuffer() :m_nCapacity(0), m_nLen(0), m_pBuf(NULL), m_pStart(NULL)
{
	init();
}
CCircBuffer::CCircBuffer(uint32_t nCapacity) :m_nLen(0), m_pBuf(NULL), m_pStart(NULL)
{
	m_nCapacity = nCapacity;
	init();
}
CCircBuffer::~CCircBuffer()
{
	if (m_pBuf != NULL)
		delete m_pBuf;
#ifdef SHOWLOG
	delete m_locallog;
#endif
}

void CCircBuffer::init()
{
	m_pBuf = new int16_t[m_nCapacity];
	m_pStart = m_pBuf;
	m_nLen = 0;
	//mutex = PTHREAD_MUTEX_INITIALIZER;
#ifdef SHOWLOG
	m_locallog = new CLocalLog("CirBuf.log");
#endif
}

int16_t CCircBuffer::CircBufReset()
{
	//pthread_mutex_lock(&mutex);
	m_pStart = m_pBuf;
	m_nLen = 0;
	//pthread_mutex_unlock(&mutex);
	return 0;
}
uint32_t CCircBuffer::CircBufGetLen()
{
	return m_nLen;
}
uint32_t CCircBuffer::CircBufGetSpace()
{
	return m_nCapacity - m_nLen;
}

int16_t CCircBuffer::CircBufSetLen(uint32_t nLen)
{
	//pthread_mutex_lock(&mutex);
	if (nLen > m_nCapacity)
	{
		//pthread_mutex_unlock(&mutex);
		return -1;
	}
	m_nLen = nLen;
	//pthread_mutex_unlock(&mutex);
	return 0;
}

int16_t CCircBuffer::CircBufAdvReadPtr(uint32_t nCount)
{
	if (nCount >= m_nLen)
		return CircBufReset();

	m_pStart += nCount;
	if (m_pStart >= (m_pBuf + m_nCapacity))
		m_pStart -= m_nCapacity;
	m_nLen -= nCount;

#ifdef SHOWLOG
	m_locallog->log(CLocalLog::T_DEBUG, "CircBuf: %d/%d", m_nLen,m_nCapacity);
#endif
	//printf("CircBuf: %d/%d \n", m_nLen, m_nCapacity);

	return 0;
}
int16_t CCircBuffer::CircBufAdvWritePtr(uint32_t nCount)
{
	if ((nCount + m_nLen) > m_nCapacity)
		return -1;

	m_nLen += nCount;

	return 0;
}

int16_t CCircBuffer::CircBufGetReadRegions(int16_t **reg1, uint32_t *reg1_len, int16_t **reg2, uint32_t *reg2_len)
{
	*reg1 = m_pStart;
	*reg1_len = m_nLen;

	if (*reg1 + *reg1_len > m_pBuf + m_nCapacity) {
		*reg1_len = (unsigned)(m_pBuf + m_nCapacity - m_pStart);
		*reg2 = m_pBuf;
		*reg2_len = m_nLen - *reg1_len;
	}
	else {
		*reg2 = NULL;
		*reg2_len = 0;
	}
	if (*reg1_len == 0)
		return -1;
	if (*reg1_len == 0 && m_nLen != 0)
		return -1;
	if ((*reg1_len + *reg2_len) != m_nLen)
		return -1;
	return 0;
}
int16_t CCircBuffer::CircBufGetWriteRegions(int16_t **reg1, uint32_t *reg1_len, int16_t **reg2, uint32_t *reg2_len)
{
	*reg1 = m_pStart + m_nLen;

	if (*reg1 >= m_pBuf + m_nCapacity)
		*reg1 -= m_nCapacity;
	*reg1_len = m_nCapacity - m_nLen;

	if (*reg1 + *reg1_len > m_pBuf + m_nCapacity) {
		*reg1_len = (unsigned)(m_pBuf + m_nCapacity - *reg1);
		*reg2 = m_pBuf;
		*reg2_len = (unsigned)(m_pStart - m_pBuf);
	}
	else {
		*reg2 = NULL;
		*reg2_len = 0;
	}
	return 0;
}

/*
count: in samples/word 16
*/
int16_t CCircBuffer::CircBufRead(int16_t *data, uint32_t count)
{
	int16_t *reg1, *reg2;
	uint32_t reg1cnt = 0;
	uint32_t reg2cnt = 0;
	int16_t	 nRet = 0;

	//pthread_mutex_lock(&mutex);

	/* Data in the buffer is less than requested */
	if (count > m_nLen)
	{
		//pthread_mutex_unlock(&mutex);
		return -1;
	}

	CircBufGetReadRegions(&reg1, &reg1cnt,&reg2, &reg2cnt);
	if (reg1cnt >= count) {
		memcpy_int16(data, reg1, count);
	}
	else {
		memcpy_int16(data, reg1, reg1cnt);
		memcpy_int16(data + reg1cnt, reg2, count - reg1cnt);
	}
	nRet = CircBufAdvReadPtr(count);

	//pthread_mutex_unlock(&mutex);

	return nRet;
}

int16_t CCircBuffer::CircBufWrite(int16_t *data, uint32_t count)
{
	int16_t *reg1, *reg2;
	uint32_t reg1cnt = 0;
	uint32_t reg2cnt = 0;
	int16_t	nRet = 0;

	//pthread_mutex_lock(&mutex);

	/* Data to write is larger than buffer can store */
	if (count > m_nCapacity - m_nLen)
	{
		//pthread_mutex_unlock(&mutex);
		return -1;
	}

	CircBufGetWriteRegions(&reg1, &reg1cnt,&reg2, &reg2cnt);
	if (reg1cnt >= count) {
		memcpy_int16(reg1, data, count);
	}
	else {
		memcpy_int16(reg1, data, reg1cnt);
		memcpy_int16(reg2, data + reg1cnt, count - reg1cnt);
	}
	nRet = CircBufAdvWritePtr(count);

	//pthread_mutex_unlock(&mutex);

	return nRet;
}

int16_t CCircBuffer::CircBufPeek(int16_t *data, uint32_t count)
{
	int16_t *reg1, *reg2;
	uint32_t reg1cnt = 0;
	uint32_t reg2cnt = 0;

	//pthread_mutex_lock(&mutex);
	/* Data in the buffer is less than requested */
	if (count > m_nLen)
	{
		//pthread_mutex_unlock(&mutex);
		return -1;
	}

	CircBufGetReadRegions(&reg1, &reg1cnt, &reg2, &reg2cnt);
	if (reg1cnt >= count) {
		memcpy_int16(data, reg1, count);
	}
	else {
		memcpy_int16(data, reg1, reg1cnt);
		memcpy_int16(data + reg1cnt, reg2, count - reg1cnt);
	}
	//pthread_mutex_unlock(&mutex);

	return 0;
}