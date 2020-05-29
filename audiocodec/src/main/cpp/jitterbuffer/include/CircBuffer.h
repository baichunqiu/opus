#ifndef _CIRC_BUFFER_H
#define _CIRC_BUFFER_H

#include "datatypes.h"
#ifdef SHOWLOG
#include "locallog.h"
#endif
#include <pthread.h>

class CCircBuffer
{
public:
	CCircBuffer();
	CCircBuffer(uint32_t nCapacity);
	~CCircBuffer();
public:
	int16_t		CircBufReset();

	uint32_t	CircBufGetLen();	//Get the data length
	uint32_t	CircBufGetSpace();	//Get the Space

	int16_t		CircBufSetLen(uint32_t nLen);

	int16_t		CircBufRead(int16_t *data, uint32_t count);
	int16_t		CircBufWrite(int16_t *data, uint32_t count);

	int16_t		CircBufPeek(int16_t *data, uint32_t count);

	int16_t		CircBufAdvReadPtr(uint32_t nCount);
	int16_t		CircBufAdvWritePtr(uint32_t nCount);

private:
	void	init();

	
	/*region 1,region 2,deal with the round trip*/
	int16_t	CircBufGetReadRegions(int16_t **reg1, uint32_t *reg1_len, int16_t **reg2, uint32_t *reg2_len);
	int16_t	CircBufGetWriteRegions(int16_t **reg1, uint32_t *reg1_len, int16_t **reg2, uint32_t *reg2_len);


private:
	int16_t		*m_pBuf;	    /**< The buffer						*/
	uint32_t	m_nCapacity;	/**< Buffer capacity, in samples	*/

	int16_t		*m_pStart;	    /**< Pointer to the first sample	*/
	uint32_t	m_nLen;		/**< Audio samples length,in samples*/

private:
#ifdef SHOWLOG
	CLocalLog	*m_locallog;
#endif
	//pthread_mutex_t mutex;
};

#endif

