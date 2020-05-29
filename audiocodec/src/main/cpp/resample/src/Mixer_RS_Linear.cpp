//====================================================================
// File Name : Mixer_RS_Linear.cpp
// Function  : Audio Mixer Resample Linear Definition file
// Program   : Li, Shaoen (LSN)
// Date      : Jun 30, 2008
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================
//#include "Mixer_Inc.h"
#include "Mixer_Def.h"

#include "Mixer_Resample.h"
#include "Mixer_RS_Linear.h"

#if	MIXER_API_EN

//	Constructor
Mixer_RS_Linear::Mixer_RS_Linear()
{
	m_nChannels			= 2;
	m_nDTB				= (1<<NP);
	m_dblRatio			= 1.0;
	m_nLastSample[0]	= 0;
	m_nLastSample[1]	= 0;
	m_nLastRatio		= 0;
	m_nNextInx			= 0;
}

//	Destructor
Mixer_RS_Linear::~Mixer_RS_Linear()
{
}

//====================================================================
// Function			: BOOL	Mixer_RS_Linear::SetMethod(UINT nMethod, UINT nChannels)
// Parameter		: UINT nMethod		- the method, defined in samplerate.h
//					  UINT nChannels	- channel number
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Method.
//====================================================================
BOOL	Mixer_RS_Linear::SetMethod(UINT, UINT nChannels)
{
	m_nChannels	= nChannels;
	return TRUE;
}

//====================================================================
// Function			: BOOL	Mixer_RS_Linear::SetRatio(UINT nOutSampleRate, UINT nInSampleRate)
// Parameter		: UINT nOutSampleRate	- OUT samplerate
//					  UINT nInSampleRate	- IN samplerate
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Ratio.
//====================================================================
BOOL	Mixer_RS_Linear::SetRatio(UINT nOutSampleRate, UINT nInSampleRate)
{
	return SetRatio((DOUBLE)nOutSampleRate / nInSampleRate);
}

//====================================================================
// Function			: BOOL	Mixer_RS_Linear::SetRatio(DOUBLE dblRatio)
// Parameter		: DOUBLE dblRatio	- OUT samplerate / IN samplerate
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Ratio.
//====================================================================
BOOL	Mixer_RS_Linear::SetRatio(DOUBLE dblRatio)
{
	m_dblRatio	= dblRatio;
	m_nDTB		= (INT)((1/m_dblRatio) * (1<<NP) + 0.5);
	return TRUE;
}

//====================================================================
// Function			: BOOL	Mixer_RS_Linear::Reset()
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Reset
//====================================================================
BOOL	Mixer_RS_Linear::Reset()
{
	m_nChannels			= 2;
	m_nDTB				= (1<<NP);
	m_dblRatio			= 1.0;
	m_nLastSample[0]	= 0;
	m_nLastSample[1]	= 0;
	m_nLastRatio		= 0;
	m_nNextInx			= 0;

	return TRUE;
}

//====================================================================
// Function			: INT		Mixer_RS_Linear::ConvertUp(SHORT* lpDes,
//					  INT nDesFNum, FLOAT* lpSrc, INT nSrcFNum, INT* pUsedNum)
// Parameter		: SHORT* lpDes		- the des buffer
//					  INT nDesFrameNum	- the des buffer frame number
//					  SHORT* lpSrc		- the src buffer
//					  INT nSrcFrameNum	- the src buffer frame number
//					  INT* pUsedFrameNum	- the src buffer used frame number
// Returned Values	: If success, return sample number; or return -1
// Memo				: Convert Up.
//====================================================================
INT		Mixer_RS_Linear::ConvertUp(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum)
{
	INT	nUsedNum	= 0; //src used
	INT	nGenNum		= 0; //des used

	INT	nVal, nX1, nX2;
	SHORT*	pSrc	= lpSrc;
	SHORT*	pDes	= lpDes;
	INT	nSrcFN		= nSrcFrameNum - 1; //src num

	if(m_nNextInx)
	{
		while(m_nNextInx < 0)
		{
			for(INT i=0; i<m_nChannels; i++)
			{
				nX1			= m_nLastSample[i];
				nX2			= pSrc[i];

				nVal		= nX1 * ((1<<NP) - m_nLastRatio) + nX2 * m_nLastRatio;
				pDes[i]		= IntToShort(nVal);
			}

			pDes	+= m_nChannels;
			nGenNum++;

			m_nLastRatio	+= m_nDTB;
			m_nNextInx	+= (m_nLastRatio >> NP);
			m_nLastRatio	&= PMASK;
		}

		pSrc	+= m_nChannels * m_nNextInx;
		nUsedNum+= m_nNextInx;
	}

	while(TRUE)
	{
		for(INT i=0; i<m_nChannels; i++)
		{
			nX1			= pSrc[i];
			nX2			= pSrc[i+m_nChannels];

			nVal		= nX1 * ((1<<NP) - m_nLastRatio) + nX2 * m_nLastRatio;
			pDes[i]		= IntToShort(nVal);
		}

		pDes	+= m_nChannels;
		nGenNum++;

		m_nLastRatio	+= m_nDTB;

		pSrc	+= m_nChannels * (m_nLastRatio >> NP);
		nUsedNum+= (m_nLastRatio >> NP);

		if(nUsedNum >= nSrcFN)
		{
			if(nUsedNum == nSrcFN)
			{
				for(INT i=0; i<m_nChannels; i++)
					m_nLastSample[i]	= pSrc[i];

				m_nNextInx	= -1;
			}
			else
			{
				m_nNextInx	= nUsedNum - nSrcFN - 1;
			}

			m_nLastRatio	&= PMASK;
			break;
		}
		m_nLastRatio	&= PMASK;
	}

	ASSERT(nGenNum <= nDesFrameNum);
	if(pUsedFrameNum)
		*pUsedFrameNum	= nUsedNum + 1;
	return nGenNum;
}

//====================================================================
// Function			: INT		Mixer_RS_Linear::ConvertDown(SHORT* lpDes,
//					  INT nDesFNum, FLOAT* lpSrc, INT nSrcFNum, INT* pUsedNum)
// Parameter		: SHORT* lpDes		- the des buffer
//					  INT nDesFrameNum	- the des buffer frame number
//					  SHORT* lpSrc		- the src buffer
//					  INT nSrcFrameNum	- the src buffer frame number
//					  INT* pUsedFrameNum	- the src buffer used frame number
// Returned Values	: If success, return sample number; or return -1
// Memo				: Convert Down.
//====================================================================
INT		Mixer_RS_Linear::ConvertDown(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum)
{
	INT	nUsedNum	= 0;
	INT	nGenNum		= 0;

	INT	nVal, nX1, nX2;
	SHORT*	pSrc	= lpSrc;
	SHORT*	pDes	= lpDes;
	INT	nSrcFN		= nSrcFrameNum - 1;

	if(m_nNextInx)
	{
		if(m_nNextInx < 0)
		{
			for(INT i=0; i<m_nChannels; i++)
			{
				nX1			= m_nLastSample[i];
				nX2			= pSrc[i];

				nVal		= nX1 * ((1<<NP) - m_nLastRatio) + nX2 * m_nLastRatio;
				pDes[i]		= IntToShort(nVal);
			}

			pDes	+= m_nChannels;
			nGenNum++;

			m_nLastRatio	+= m_nDTB;
			m_nNextInx	+= (m_nLastRatio >> NP);
			m_nLastRatio	&= PMASK;
		}

		pSrc	+= m_nChannels * m_nNextInx;
		nUsedNum+= m_nNextInx;
	}

	while(TRUE)
	{
		for(INT i=0; i<m_nChannels; i++)
		{
			nX1			= pSrc[i];
			nX2			= pSrc[i+m_nChannels];

			nVal		= nX1 * ((1<<NP) - m_nLastRatio) + nX2 * m_nLastRatio;
			pDes[i]		= IntToShort(nVal);
		}

		pDes	+= m_nChannels;
		nGenNum++;

		m_nLastRatio	+= m_nDTB;

		pSrc	+= m_nChannels * (m_nLastRatio >> NP);
		nUsedNum+= (m_nLastRatio >> NP);

		if(nUsedNum >= nSrcFN)
		{
			if(nUsedNum == nSrcFN)
			{
				for(INT i=0; i<m_nChannels; i++)
					m_nLastSample[i]	= pSrc[i];

				m_nNextInx	= -1;
			}
			else
			{
				m_nNextInx	= nUsedNum - nSrcFN - 1;
			}

			m_nLastRatio	&= PMASK;
			break;
		}
		m_nLastRatio	&= PMASK;
	}

	ASSERT(nGenNum <= nDesFrameNum);
	if(pUsedFrameNum)
		*pUsedFrameNum	= nUsedNum + 1;
	return nGenNum;
}


//====================================================================
//	Protected Functions
//====================================================================

//====================================================================
// Function			: SHORT	Mixer_RS_Linear::IntToShort(INT nVal)
// Parameter		: INT nVal	- int value
// Returned Values	: Return the short value
// Memo				: Int to Short.
//====================================================================
SHORT	Mixer_RS_Linear::IntToShort(INT nVal)
{
	nVal	+= (1 << (NP - 1));
	nVal	>>= NP;

	if(nVal > MAX_SHORT)
		return (SHORT)MAX_SHORT;
	if(nVal < MIN_SHORT)
		return (SHORT)MIN_SHORT;

	return (SHORT)nVal;
}

#endif
