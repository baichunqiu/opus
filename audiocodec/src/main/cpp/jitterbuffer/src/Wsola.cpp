#include "Wsola.h"
#include <math.h>

CWsola::CWsola() :m_nSamplingRate(48000), m_nChannel(2), m_nFramesInterval(20), m_nWindow_Size(64)
{
	m_nSamplesPerFrame = m_nSamplingRate * m_nFramesInterval * m_nChannel / 1000;

	m_nTempBufSize = m_nSamplesPerFrame * 10;/*we use 10 frames to buffer the input data*/
	m_nExpand_Max_Dist = 1.5*m_nSamplesPerFrame / m_nChannel;
	m_nExpand_Min_Dist = 0.5*m_nSamplesPerFrame / m_nChannel;
}
CWsola::~CWsola()
{
	if (m_pHanning_Window != NULL)
		delete m_pHanning_Window;
	if (m_pMergeBuf != NULL)
		delete m_pMergeBuf;
	if (m_pTempExpandL != NULL)
		delete m_pTempExpandL;
	if (m_pTempExpandR != NULL)
		delete m_pTempExpandR;
	if (m_pTempCompressL != NULL)
		delete m_pTempCompressL;
	if (m_pTempCompressR != NULL)
		delete m_pTempCompressR;

#ifdef SHOWLOG
	if (m_locallog != NULL)
		delete m_locallog;
#endif
}

CWsola::CWsola(uint16_t nSamplingRate, uint16_t nChannel,
	uint16_t nFramesInterval,
	uint16_t nWindowPitchTime)
{
	uint16_t nWindow_Size = 0;
	m_nSamplingRate = nSamplingRate;
	m_nChannel = nChannel;
	m_nFramesInterval = nFramesInterval;

	nWindow_Size = m_nSamplingRate * nWindowPitchTime / 1000;
	m_nSamplesPerFrame = m_nSamplingRate * m_nFramesInterval * m_nChannel / 1000;

	m_nWindow_Size = nWindow_Size;//Hanning Window Size


	m_nTemplateSize = nWindow_Size;
	if (m_nTemplateSize > m_nSamplesPerFrame)
		m_nTemplateSize = m_nSamplesPerFrame;

	m_nTempBufSize = m_nSamplesPerFrame * 10;/*we use 10 frames to buffer the input data*/

	m_nExpand_Max_Dist = 1.5*m_nSamplesPerFrame / m_nChannel; //every channel we need to search
	m_nExpand_Min_Dist = 0.5*m_nSamplesPerFrame / m_nChannel; //every channel we need to search

	init();
}
void CWsola::init()
{
	m_pHanning_Window = new float[m_nWindow_Size];
	m_pMergeBuf = new int16_t[m_nWindow_Size];

	m_pTempExpandL = new int16_t[m_nTempBufSize];//this is used for data tempory putted in
	m_pTempExpandR = new int16_t[m_nTempBufSize];

	m_pTempCompressL = new int16_t[m_nTempBufSize];
	m_pTempCompressR = new int16_t[m_nTempBufSize];

	create_win(m_pHanning_Window, m_nWindow_Size);

#ifdef SHOWLOG
	m_locallog = new CLocalLog("WSola.log");
#endif
}
/*
buf: the int16 pcm buffer,original buffer,this buffer at least contain 2 frames of data
and the buffer length at least 3 frames of data
nLen: in int16_t words*/

/*nExpandedCnt = m_pWsola->wsola_generate(m_pExpandBuf, nPCMLen, m_nExpandBufSize, m_nSamplesPerFrame);
*/
int16_t CWsola::wsola_generate(int16_t buf[], uint32_t nPCMLen, uint32_t nBufLen, uint32_t nExpandedLen_ST)
{
	//Do some sanity check
	uint16_t nExpandedLen = 0;
	nExpandedLen = nExpandedLen_ST / m_nChannel;
	if (nExpandedLen > m_nExpand_Max_Dist)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: too more data to be generated");
#endif
		return -1;
	}
	if (nBufLen <= (nPCMLen + nExpandedLen_ST))
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: buf is not enough,expand can not be executed");
#endif
		return -1;
	}
	if ((nPCMLen + nExpandedLen_ST) > m_nTempBufSize * m_nChannel)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: temp buf size is too small");
#endif
		return -1;
	}
	if (m_nTempBufSize < (nPCMLen / m_nChannel + m_nExpand_Max_Dist))
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: temp buf size is not enough");
#endif
		return -1;
	}

	uint16_t nChannelLen = 0;
	uint16_t nIndex = 0;
	uint16_t nGenLen = 0;

	if (m_nChannel == 2)
	{
		nChannelLen = nPCMLen / 2;
		for (nIndex = 0; nIndex < nChannelLen; nIndex++)
		{
			m_pTempExpandL[nIndex] = buf[2 * nIndex];
			m_pTempExpandR[nIndex] = buf[2 * nIndex + 1];
		}
		/*We should perform two expand operation*/
		nGenLen = expand(m_pTempExpandL, m_pTempExpandR, nChannelLen, nExpandedLen);
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "wsola_generate:: %d/%d", nGenLen, nExpandedLen);
#endif
		/*copy the expanded data into original buffer*/
		for (nIndex = 0; nIndex < nChannelLen + nGenLen; nIndex++)
		{
			buf[2 * nIndex] = m_pTempExpandL[nIndex];
			buf[2 * nIndex + 1] = m_pTempExpandR[nIndex];
		}
		return nGenLen * 2;

	}
	else if (m_nChannel == 1)
	{
		nChannelLen = nPCMLen;
		for (nIndex = 0; nIndex < nChannelLen; nIndex++)
		{
			m_pTempExpandL[nIndex] = buf[nIndex];
		}
		nGenLen = expand(m_pTempExpandL, nChannelLen, nExpandedLen);
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "wsola_generate:: %d", nGenLen);
#endif
		for (nIndex = 0; nIndex < nChannelLen + nGenLen; nIndex++)
		{
			buf[nIndex] = m_pTempExpandL[nIndex];
		}
		return nGenLen;
	}
	else
		return -1;

	return -1;
}
int16_t CWsola::wsola_discard(int16_t buf[], uint32_t nPCMLen, uint32_t nBufLen, uint32_t nCompressedLen_ST)
{
	//Do some sanity check
	if (nPCMLen < (nCompressedLen_ST + m_nWindow_Size * m_nChannel))
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: buf data is not enough");
#endif
		return -1;
	}
	if (nPCMLen > nBufLen)
	{
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "Wsola:: this can not be true");
#endif
		return -1;
	}

	uint32_t nChannelLen = 0;
	uint16_t nIndex = 0;
	uint16_t nComLenL = 0;
	uint16_t nComLenR = 0;
	uint16_t nComLen = 0;
	uint32_t nCompressedLen = 0;

	uint32_t nForCompressToReadLen = 0;

	nCompressedLen = nCompressedLen_ST / m_nChannel;

	if (m_nChannel == 2)
	{
		nChannelLen = nPCMLen / 2;
		//nChannelLen: single channel data length
		//m_nTempBufSize: single channel data buf size
		if (nChannelLen >= m_nTempBufSize)
		{
			nForCompressToReadLen = m_nTempBufSize;
		}
		else
		{
			nForCompressToReadLen = nChannelLen;
		}

		for (nIndex = 0; nIndex < nForCompressToReadLen; nIndex++)
		{
			m_pTempCompressL[nIndex] = buf[2 * nIndex];
			m_pTempCompressR[nIndex] = buf[2 * nIndex + 1];
		}
		/*We should perform two expand operation*/
		/*for compress,we compress the begining*/
		nComLen = compress(m_pTempCompressL, m_pTempCompressR, nForCompressToReadLen, nCompressedLen);
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "wsola_discard:: %d/%d", nComLen, nCompressedLen);
#endif
		if (nForCompressToReadLen <= nChannelLen)
		{
			//perform a memory move and write
			memmove_int16(buf + 2 * (nForCompressToReadLen - nComLen), buf + 2 * nForCompressToReadLen, 
				2 * (nChannelLen - nForCompressToReadLen));

			for (nIndex = 0; nIndex < nForCompressToReadLen - nComLen; nIndex++)
			{
				buf[2 * nIndex] = m_pTempCompressL[nIndex];
				buf[2 * nIndex + 1] = m_pTempCompressR[nIndex];
			}
		}
		else
		{
			/*copy the expanded data into original buffer*/
			for (nIndex = 0; nIndex < nForCompressToReadLen - nComLen; nIndex++)
			{
				buf[2 * nIndex] = m_pTempCompressL[nIndex];
				buf[2 * nIndex + 1] = m_pTempCompressR[nIndex];
			}
		}
		
		return nComLen * 2;

	}
	else if (m_nChannel == 1)
	{
		nChannelLen = nPCMLen;
		for (nIndex = 0; nIndex < nChannelLen; nIndex++)
		{
			m_pTempCompressL[nIndex] = buf[nIndex];
		}
		nComLen = compress(m_pTempCompressL, nChannelLen, nCompressedLen);
#ifdef SHOWLOG
		m_locallog->log(CLocalLog::T_DEBUG, "wsola_discard:: %d", nComLen);
#endif
		for (nIndex = 0; nIndex < nChannelLen - nComLen; nIndex++)
		{
			buf[nIndex] = m_pTempCompressL[nIndex];
		}
		return nComLen;
	}
	else
		return -1;

	return -1;
}
/**/
uint16_t CWsola::compress(int16_t *buf, uint16_t count, uint16_t del_cnt)
{
	unsigned samples_del = 0, rep;

	for (rep = 1;; ++rep) {
		int16_t *start, *end;
		unsigned dist;

		if (count <= m_nWindow_Size + del_cnt) {
			return samples_del;
		}

		// Make start distance to del_cnt, so discard will be performed in
		// only one iteration.
		//start = buf + (frmsz >> 1);
		start = buf + del_cnt - samples_del;
		end = start + m_nSamplesPerFrame;

		if (end + m_nWindow_Size > buf + count) {
			end = buf + count - m_nWindow_Size;
		}

		start = find_pitch(buf, start, end, m_nTemplateSize, 0);
		dist = (unsigned)(start - buf);

		overlapp_add(buf, m_nWindow_Size, buf, start, m_pHanning_Window);

		memmove_int16(buf + m_nWindow_Size,
			buf + m_nWindow_Size + dist,
			count - m_nWindow_Size - dist);

		count -= dist;
		samples_del += dist;

		if (samples_del >= del_cnt) {
			break;
		}
	}

	return samples_del;
}
uint16_t CWsola::compress(int16_t *pcm_buf_l, int16_t *pcm_buf_r, uint16_t pcm_size, uint16_t del_cnt)
{
	unsigned samples_del = 0, rep;
	int16_t *startl, *endl;
	int16_t *startr, *endr;
	double pitch_corr_val_l, pitch_corr_val_r;

	for (rep = 1;; ++rep) {

		unsigned dist;

		if (pcm_size <= m_nWindow_Size + del_cnt) {
			return samples_del;
		}

		// Make start distance to del_cnt, so discard will be performed in
		// only one iteration.
		//start = buf + (frmsz >> 1);
		//////////////////////////////////////////////////////////////////////////////
		startl = pcm_buf_l + del_cnt - samples_del;
		endl = startl + m_nSamplesPerFrame;

		if (endl + m_nWindow_Size > pcm_buf_l + pcm_size) {
			endl = pcm_buf_l + pcm_size - m_nWindow_Size;
		}
		startr = pcm_buf_r + del_cnt - samples_del;
		endr = startr + m_nSamplesPerFrame;

		if (endr + m_nWindow_Size > pcm_buf_r + pcm_size) {
			endr = pcm_buf_r + pcm_size - m_nWindow_Size;
		}
		//////////////////////////////////////////////////////////////////////////////
		startl = find_pitch(pcm_buf_l, startl, endl, m_nTemplateSize, 0, &pitch_corr_val_l);
		startr = find_pitch(pcm_buf_r, startr, endr, m_nTemplateSize, 0, &pitch_corr_val_r);
		if (pitch_corr_val_l > pitch_corr_val_r)
		{
			startl = startl;
			startr = (unsigned)(startl - pcm_buf_l) + pcm_buf_r;
		}
		else
		{
			startl = (unsigned)(startr - pcm_buf_r) + pcm_buf_l;
			startr = startr;
		}
		//////////////////////////////////////////////////////////////////////////////
		dist = (unsigned)(startl - pcm_buf_l);

		overlapp_add(pcm_buf_l, m_nWindow_Size, pcm_buf_l, startl, m_pHanning_Window);

		memmove_int16(pcm_buf_l + m_nWindow_Size,
			pcm_buf_l + m_nWindow_Size + dist,
			pcm_size - m_nWindow_Size - dist);
		//////////////////////////////////////////////////////////////////////////////
		dist = (unsigned)(startr - pcm_buf_r);

		overlapp_add(pcm_buf_r, m_nWindow_Size, pcm_buf_r, startr, m_pHanning_Window);

		memmove_int16(pcm_buf_r + m_nWindow_Size,
			pcm_buf_r + m_nWindow_Size + dist,
			pcm_size - m_nWindow_Size - dist);
		//////////////////////////////////////////////////////////////////////////////

		pcm_size -= dist;
		samples_del += dist;

		if (samples_del >= del_cnt) {
			break;
		}
	}

	return samples_del;
}

/*
	nGenLen = expand(m_pTempExpandL, m_pTempExpandR, nChannelLen, nExpandedLen);
*/
uint16_t CWsola::expand(int16_t *pcm_buf_l, int16_t *pcm_buf_r, uint16_t pcm_size, uint16_t needed)
{
	uint16_t generated = 0;
	uint16_t rep;
	int16_t *startl, *startr, *templ, *tempr;
	double pitch_corr_val_l, pitch_corr_val_r;

	for (rep = 1;; ++rep) {

		unsigned dist;

		templ = pcm_buf_l + pcm_size - m_nWindow_Size;
		tempr = pcm_buf_r + pcm_size - m_nWindow_Size;

		startl = find_pitch(templ,
			templ - m_nExpand_Max_Dist,
			templ - m_nExpand_Min_Dist,
			m_nTemplateSize,
			0, //the first one or the best one
			&pitch_corr_val_l);

		startr = find_pitch(tempr,
			tempr - m_nExpand_Max_Dist,
			tempr - m_nExpand_Min_Dist,
			m_nTemplateSize,
			0, //the first one or the best one
			&pitch_corr_val_r);
		if (pitch_corr_val_l > pitch_corr_val_r)
		{
			//We will use the left speficied position
			startl = startl;
			startr = tempr - (unsigned)(templ - startl);
		}
		else
		{
			//We will use the left speficied position
			startl = templ - (unsigned)(tempr - startr);
			startr = startr;
		}
		////////////////////////////////////////////////
		overlapp_add(m_pMergeBuf, m_nWindow_Size, templ,
			startl, m_pHanning_Window);
		/* How many new samples do we have */
		dist = (unsigned)(templ - startl);
		/* Copy the "tail" (excess frame) to the end */
		memmove_int16(templ + m_nWindow_Size,
			startl + m_nWindow_Size,
			dist);
		/* Copy the merged frame */
		memcpy_int16(templ, m_pMergeBuf, m_nWindow_Size);
		////////////////////////////////////////////////
		overlapp_add(m_pMergeBuf, m_nWindow_Size, tempr,
			startr, m_pHanning_Window);
		/* How many new samples do we have */
		dist = (unsigned)(tempr - startr);
		/* Copy the "tail" (excess frame) to the end */
		memmove_int16(tempr + m_nWindow_Size,
			startr + m_nWindow_Size,
			dist);
		/* Copy the merged frame */
		memcpy_int16(tempr, m_pMergeBuf, m_nWindow_Size);
		////////////////////////////////////////////////
		generated += dist;
		pcm_size += dist;

		if (generated >= needed) {
			break;
		}
	}
	return generated;
}

uint16_t CWsola::expand(int16_t *pcm_buf, uint16_t buf_size, uint16_t needed)
{
	uint16_t generated = 0;
	uint16_t rep;

	for (rep = 1;; ++rep) {
		int16_t *start, *templ;
		unsigned dist;

		templ = pcm_buf + buf_size - m_nWindow_Size;

		start = find_pitch(templ,
			templ - m_nExpand_Max_Dist,
			templ - m_nExpand_Min_Dist,
			m_nTemplateSize,
			0 //the first one or the best one
			);

		overlapp_add(m_pMergeBuf, m_nWindow_Size, templ,
			start, m_pHanning_Window);

		/* How many new samples do we have */
		dist = (unsigned)(templ - start);

		/* Copy the "tail" (excess frame) to the end */
		memmove_int16(templ + m_nWindow_Size,
			start + m_nWindow_Size,
			dist);

		/* Copy the merged frame */
		memcpy_int16(templ, m_pMergeBuf, m_nWindow_Size);

		generated += dist;

		buf_size += dist;

		if (generated >= needed) {
			break;
		}
	}
	return generated;
}

int16_t *CWsola::find_pitch(int16_t *from, int16_t *begin, int16_t *end,
	uint16_t template_size, int16_t first)
{
	int16_t *sr, *best = begin;
	double best_corr = 0;

	for (sr = begin; sr != end; ++sr) {
		double corr = 0;
		uint16_t i;

		/* Do calculation on 8 samples at once */
		for (i = 0; i < template_size - 8; i += 8) {
			corr += ((float)from[i + 0]) * ((float)sr[i + 0]) +
				((float)from[i + 1]) * ((float)sr[i + 1]) +
				((float)from[i + 2]) * ((float)sr[i + 2]) +
				((float)from[i + 3]) * ((float)sr[i + 3]) +
				((float)from[i + 4]) * ((float)sr[i + 4]) +
				((float)from[i + 5]) * ((float)sr[i + 5]) +
				((float)from[i + 6]) * ((float)sr[i + 6]) +
				((float)from[i + 7]) * ((float)sr[i + 7]);
		}

		/* Process remaining samples. */
		for (; i < template_size; ++i) {
			corr += ((float)from[i]) * ((float)sr[i]);
		}

		if (first) {
			if (corr > best_corr) {
				best_corr = corr;
				best = sr;
			}
		}
		else {
			if (corr >= best_corr) {
				best_corr = corr;
				best = sr;
			}
		}
	}

	return best;
}

/*startl = find_pitch(templ,
	templ - m_nExpand_Max_Dist,
	templ - m_nExpand_Min_Dist,
	m_nTemplateSize,
	0, //the first one or the best one
	&pitch_corr_val_l);
*/

int16_t *CWsola::find_pitch(int16_t *from, int16_t *begin, int16_t *end,
	uint16_t template_size, int16_t first, double* corr_value)
{
	int16_t *sr, *best = begin;
	double best_corr = 0;

	for (sr = begin; sr != end; ++sr) {
		double corr = 0;
		uint16_t i;

		/* Do calculation on 8 samples at once */
		for (i = 0; i < template_size - 8; i += 8) {
			corr += ((float)from[i + 0]) * ((float)sr[i + 0]) +
				((float)from[i + 1]) * ((float)sr[i + 1]) +
				((float)from[i + 2]) * ((float)sr[i + 2]) +
				((float)from[i + 3]) * ((float)sr[i + 3]) +
				((float)from[i + 4]) * ((float)sr[i + 4]) +
				((float)from[i + 5]) * ((float)sr[i + 5]) +
				((float)from[i + 6]) * ((float)sr[i + 6]) +
				((float)from[i + 7]) * ((float)sr[i + 7]);
		}

		/* Process remaining samples. */
		for (; i < template_size; ++i) {
			corr += ((float)from[i]) * ((float)sr[i]);
		}

		if (first) {
			if (corr > best_corr) {
				best_corr = corr;
				best = sr;
				*corr_value = best_corr;
			}
		}
		else {
			if (corr >= best_corr) {
				best_corr = corr;
				best = sr;
				*corr_value = best_corr;
			}
		}
	}

	return best;
}

void CWsola::overlapp_add(int16_t dst[], uint16_t window_size,
	int16_t templ[], int16_t signal_match[], float window[])
{
	unsigned i;

	for (i = 0; i < window_size; ++i) {
		dst[i] = (int16_t)(templ[i] * window[window_size - 1 - i] + signal_match[i] * window[i]);
	}
}

void CWsola::create_win(float *pw, uint16_t window_size)
{
	uint16_t i;

	for (i = 0; i < window_size; i++) {
		pw[i] = (float)(0.5 - 0.5 * cos(2.0 * PI * i / (window_size * 2 - 1)));
	}
}