//====================================================================
// File Name : Mixer_RS_KW.cpp
// Function  : Audio Mixer Resample Kaiser-Windowed Definition file
// Program   : Li, Shaoen (LSN)
// Date      : Jun 30, 2008
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#include "Mixer_Def.h"

#include "Mixer_Resample.h"
#include "Mixer_RS_Linear.h"
#include "Mixer_RS_KW.h"

#include "largefilter.h"
#include "smallfilter.h"

#include <android/log.h>
#define TAG    "Mixer"
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG,__VA_ARGS__)

#if	1

//	Constructor
 Mixer_RS_KW::Mixer_RS_KW()
{
	m_lpIMP		= NULL;
	m_lpIMPD	= NULL;
	m_nWing		= 0;
	m_nSCL		= 0;
	m_nMult		= 0;
	m_bInterp	= FALSE;
	m_nDHB		= 0;
	m_nXOff		= 16;
}

//	Destructor
 Mixer_RS_KW::~Mixer_RS_KW()
{
}

//====================================================================
// Function			: BOOL	Mixer_RS_KW::SetMethod(UINT nMethod, UINT nChannels)
// Parameter		: UINT nMethod		- the method, defined in samplerate.h
//					  UINT nChannels	- channel number
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Method.
//====================================================================
 BOOL	Mixer_RS_KW::SetMethod(UINT nMethod, UINT nChannels)
{
	m_bInterp	= nMethod & 0xFFFF0000;
	switch(nMethod & 0xFFFF)
	{
	case 0:
		{
			m_nMult		= SMALL_FILTER_NMULT;
			m_lpIMP		= SMALL_FILTER_IMP;
			m_lpIMPD	= SMALL_FILTER_IMPD;
			m_nSCL		= SMALL_FILTER_SCALE;
			m_nWing		= SMALL_FILTER_NWING;
			break;
		}
	default:
		{
			m_nMult		= LARGE_FILTER_NMULT;
			m_lpIMP		= LARGE_FILTER_IMP;
			m_lpIMPD	= LARGE_FILTER_IMPD;
			m_nSCL		= LARGE_FILTER_SCALE;
			m_nWing		= LARGE_FILTER_NWING;
			break;
		}
	}

	m_nSCL			*= 19;
	m_nSCL			/= 20;
	return Mixer_RS_Linear::SetMethod(nMethod, nChannels);
}

//====================================================================
// Function			: BOOL	Mixer_RS_KW::SetRatio(UINT nOutSampleRate, UINT nInSampleRate)
// Parameter		: UINT nOutSampleRate	- OUT samplerate
//					  UINT nInSampleRate	- IN samplerate
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Ratio.
//====================================================================
 BOOL	Mixer_RS_KW::SetRatio(UINT nOutSampleRate, UINT nInSampleRate)
{
	return SetRatio((DOUBLE)nOutSampleRate / nInSampleRate);
}

//====================================================================
// Function			: BOOL	Mixer_RS_KW::SetRatio(DOUBLE dblRatio)
// Parameter		: DOUBLE dblRatio	- OUT samplerate / IN samplerate
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Sample Ratio.
//====================================================================
 BOOL	Mixer_RS_KW::SetRatio(DOUBLE dblRatio)
{
	m_dblRatio	= dblRatio;
	m_nDTB		= (INT)((1/m_dblRatio) * (1<<NP) + 0.5);
    m_nDHB		= (INT)(MIN(NPC, m_dblRatio*NPC) * (1<<NA) + 0.5);		/* Fixed-point representation */
	return TRUE;
}

//====================================================================
// Function			: BOOL	Mixer_RS_KW::Reset()
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Buffer Length.
//====================================================================
BOOL	Mixer_RS_KW::Reset()
{
	memset(m_aryPCM.GetDataBuf(), 0, m_nXOff * 2 * sizeof(SHORT) * m_nChannels);
	return TRUE;
}

//====================================================================
// Function			: BOOL	Mixer_RS_KW::SetPCMBufLen(UINT nBufLen)
// Parameter		: UINT nBufLen		- The local Buffer length
// Returned Values	: If success, return TRUE; or return FALSE
// Memo				: Set Buffer Length.
//====================================================================
 BOOL	Mixer_RS_KW::SetPCMBufLen(UINT nBufLen)
{
	if(!m_aryPCM.Init((nBufLen >> 1) + m_nXOff * 2 * m_nChannels))
		return FALSE;

	if(m_dblRatio < 1.0)	//For Convert Down
		m_nSCL	= (INT)(m_nSCL * m_dblRatio + 0.5);

	m_nXOff		= (INT)(((m_nMult + 1) / 2.0) * MAX(1.0, 1.0/m_dblRatio) + 10);
	memset(m_aryPCM.GetDataBuf(), 0, m_nXOff * 2 * sizeof(SHORT) * m_nChannels);
	return TRUE;
}

//====================================================================
// Function			: INT		Mixer_RS_KW::ConvertUp(SHORT* lpDes,
//					  INT nDesFNum, FLOAT* lpSrc, INT nSrcFNum, INT* pUsedNum)
// Parameter		: SHORT* lpDes		- the des buffer
//					  INT nDesFrameNum	- the des buffer frame number
//					  SHORT* lpSrc		- the src buffer
//					  INT nSrcFrameNum	- the src buffer frame number
//					  INT* pUsedFrameNum	- the src buffer used frame number
// Returned Values	: If success, return sample number; or return -1
// Memo				: Convert Up.
//====================================================================
 INT		Mixer_RS_KW::ConvertUp(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum)
{
	INT	nUsedNum	= 0;
	INT	nGenNum		= 0;

	INT	nStart		= m_nXOff * 2 * m_nChannels;

	memcpy(&m_aryPCM[nStart], lpSrc, nSrcFrameNum * m_nChannels * sizeof(SHORT));
	INT	nVal;
	SHORT*	pSrc	= &m_aryPCM[nStart/2];
	SHORT*	pDes	= lpDes;
	INT	nSrcFN		= nSrcFrameNum;

	while(TRUE)
	{

		for(INT i=0; i<m_nChannels; i++)
		{
			/* Perform left-wing inner product */
			nVal		= FilterUp(&pSrc[i], m_nLastRatio, -1);
			/* Perform right-wing inner product */
		//	nVal		+= FilterUp(&pSrc[i+m_nChannels], ((m_nLastRatio ^ PMASK) + 1) & PMASK, 1);
			nVal		+= FilterUp(&pSrc[i+m_nChannels], (1<<NP) - m_nLastRatio, 1);
			/* Make guard bits */
			nVal		>>= NHG;
			/* Normalize for unity filter gain */
			nVal		*= m_nSCL;
			pDes[i]		= IntToShort(nVal);
		}

		pDes	+= m_nChannels;
		nGenNum++;

		m_nLastRatio	+= m_nDTB;
		if(m_nLastRatio >= (1<<NP))
		{
			pSrc	+= m_nChannels;
			nUsedNum++;
			m_nLastRatio	&= PMASK;
//			avl_printf("nUsedNum = %d %d\r\n", nUsedNum, nSrcFN);
			if(nUsedNum >= nSrcFN)
			{
//			    avl_printf("nSrcFN 0x%x  m_nChannels 0x%x nStart 0x%x\r\n ", nSrcFN, m_nChannels, nStart);
				//	Reused the last values
				memcpy(&m_aryPCM[0], &m_aryPCM[nSrcFN*m_nChannels], nStart * sizeof(SHORT));
//				avl_printf("  copy done\r\n");
				break;
			}
		}
	}
//	avl_printf("output len = %d bytes\r\n", (pDes - lpDes)*2);
//	avl_printf("assert done\r\n");
	ASSERT(nGenNum <= nDesFrameNum);
	if(pUsedFrameNum)
		*pUsedFrameNum	= nUsedNum + 1;
	return nGenNum;
}

//====================================================================
// Function			: INT		Mixer_RS_KW::ConvertDown(SHORT* lpDes,
//					  INT nDesFNum, FLOAT* lpSrc, INT nSrcFNum, INT* pUsedNum)
// Parameter		: SHORT* lpDes		- the des buffer
//					  INT nDesFrameNum	- the des buffer frame number
//					  SHORT* lpSrc		- the src buffer
//					  INT nSrcFrameNum	- the src buffer frame number
//					  INT* pUsedFrameNum	- the src buffer used frame number
// Returned Values	: If success, return sample number; or return -1
// Memo				: Convert Down.
//====================================================================
 INT		Mixer_RS_KW::ConvertDown(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum)
{
	INT	nUsedNum	= 0;
	INT	nGenNum		= 0;

	INT	nStart		= m_nXOff * 2 * m_nChannels;
	memcpy(&m_aryPCM[nStart], lpSrc, nSrcFrameNum * m_nChannels * sizeof(SHORT));

	INT	nVal;
	SHORT*	pSrc	= &m_aryPCM[nStart/2];
	SHORT*	pDes	= lpDes;
	INT	nSrcFN		= nSrcFrameNum;

	while(TRUE)
	{
		for(INT i=0; i<m_nChannels; i++)
		{
			/* Perform left-wing inner product */
			nVal		= FilterDown(&pSrc[i], m_nLastRatio, -1);
			/* Perform right-wing inner product */
		//	nVal		+= FilterDown(pSrc+i+m_nChannels, (m_nLastRatio ^ PMASK) + 1, 1);
			nVal		+= FilterDown(&pSrc[i+m_nChannels], (1<<NP) - m_nLastRatio, 1);
			/* Make guard bits */
			nVal		>>= NHG;
			/* Normalize for unity filter gain */
			nVal		*= m_nSCL;
			pDes[i]		= IntToShort(nVal);
		}

		pDes	+= m_nChannels;
		nGenNum++;

		m_nLastRatio	+= m_nDTB;
		if((m_nLastRatio >= (1<<NP)) || (nUsedNum >= nSrcFN))
		{
			pSrc	+= m_nChannels * (m_nLastRatio >> NP);
			nUsedNum+= (m_nLastRatio >> NP);

			if(nUsedNum >= nSrcFN)
			{
				//	Reused the last values
				memcpy(&m_aryPCM[0], &m_aryPCM[nSrcFN * m_nChannels], nStart * sizeof(SHORT));
				break;
			}

			m_nLastRatio	&= PMASK;
		}
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
// Function			: SHORT	Mixer_RS_KW::IntToShort(INT nVal)
// Parameter		: INT nVal	- int value
// Returned Values	: Return the short value
// Memo				: Int to Short.
//====================================================================
 SHORT	Mixer_RS_KW::IntToShort(INT nVal)
{
	nVal	+= (1 << (NLPSCL - 1));
	nVal	>>= NLPSCL;

	if(nVal > MAX_SHORT)
		return (SHORT)MAX_SHORT;
	if(nVal < MIN_SHORT)
		return (SHORT)MIN_SHORT;

	return (SHORT)nVal;
}

//====================================================================
// Function			: INT	Mixer_RS_KW::FilterUp(SHORT* lpSrc,
//					  INT nPhase, INT nInc)
// Parameter		: SHORT* lpSrc		- The src data
//					  INT nPhase		- The Bias Fhase
//					  INT nInc			- The Increasemment
// Returned Values	: Return the short value
// Memo				: Int to Short.
//====================================================================
 INT		Mixer_RS_KW::FilterUp(SHORT* Xp, INT nPhase, INT nInc)
{
	SHORT	*Hp, *Hdp = NULL, *End;
	SHORT	a = 0;
	INT		v, t;

	v=0;
	Hp	= &m_lpIMP[nPhase>>NA];
	End	= &m_lpIMP[m_nWing];
	if(m_bInterp)
	{
		Hdp = &m_lpIMPD[nPhase>>NA];
		a = nPhase & AMASK;
	}
	if(nInc == 1)			/* If doing right wing...              */
	{						/* ...drop extra coeff, so when Ph is  */
		End--;				/*    0.5, we don't do too many mult's */
		if (nPhase == 0)		/* If the phase is zero...           */
		{					/* ...then we've already skipped the */
			Hp	+= NPC;		/*    first sample, so we must also  */
			Hdp	+= NPC;		/*    skip ahead in Imp[] and ImpD[] */
		}
	}
	if (m_bInterp)
	{
		while(Hp < End)
		{
			t	= *Hp;					/* Get filter coeff */
			t	+= (((INT)*Hdp)*a)>>NA;	/* t is now interp'd filter coeff */
			Hdp	+= NPC;					/* Filter coeff differences step */
			t	*= *Xp;					/* Mult coeff by input sample */
			if(t & (1<<(NHXN-1)))		/* Round, if needed */
				t	+= (1<<(NHXN-1));
			t	>>= NHXN;				/* Leave some guard bits, but come back some */
			v	+= t;					/* The filter output */
			Hp	+= NPC;					/* Filter coeff step */
			Xp	+= nInc*m_nChannels;	/* Input signal step. NO CHECK ON BOUNDS */
		}
	}
	else
	{
		while (Hp < End)
		{
			t	= *Hp;					/* Get filter coeff */
			t	*= *Xp;					/* Mult coeff by input sample */
			if(t & (1<<(NHXN-1)))		/* Round, if needed */
				t += (1<<(NHXN-1));
			t	>>= NHXN;				/* Leave some guard bits, but come back some */
			v	+= t;					/* The filter output */
			Hp	+= NPC;					/* Filter coeff step */
			Xp	+= nInc*m_nChannels;	/* Input signal step. NO CHECK ON BOUNDS */
		}
	}
	return v;
}

//====================================================================
// Function			: INT	Mixer_RS_KW::FilterDown(SHORT* lpSrc,
//					  INT nPhase, INT nInc)
// Parameter		: SHORT* lpSrc		- The src data
//					  INT nPhase		- The Bias Fhase
//					  INT nInc			- The Increasemment
// Returned Values	: Return the short value
// Memo				: Int to Short.
//====================================================================
 INT		Mixer_RS_KW::FilterDown(SHORT* Xp, INT nPhase, INT nInc)
{
	SHORT	a;
	SHORT	*Hp, *Hdp, *End;
	INT		v, t;
	UINT	Ho;

	v=0;
	Ho	= (nPhase*(UINT)m_nDHB) >> NP;
	End	= &m_lpIMP[m_nWing];

	if(nInc == 1)			/* If doing right wing...              */
	{						/* ...drop extra coeff, so when Ph is  */
		End--;				/*    0.5, we don't do too many mult's */
		if(nPhase == 0)		/* If the phase is zero...           */
			Ho += m_nDHB;	/* ...then we've already skipped the */
	}						/*    first sample, so we must also  */
							/*    skip ahead in Imp[] and ImpD[] */
	if(m_bInterp)
	{
		while((Hp = &m_lpIMP[Ho>>NA]) < End)
		{
			t	= *Hp;					/* Get IR sample */
			Hdp	= &m_lpIMPD[Ho>>NA];	/* get interp (lower Na) bits from diff table*/
			a	= Ho & AMASK;			/* a is logically between 0 and 1 */
			t	+= (((INT)*Hdp)*a)>>NA;	/* t is now interp'd filter coeff */
			t	*= *Xp;					/* Mult coeff by input sample */
			if(t & 1<<(NHXN-1))			/* Round, if needed */
				t	+= 1<<(NHXN-1);
			t	>>= NHXN;				/* Leave some guard bits, but come back some */
			v	+= t;					/* The filter output */
			Ho	+= m_nDHB;				/* IR step */
			Xp	+= nInc * m_nChannels;	/* Input signal step. NO CHECK ON BOUNDS */
		}
	}
	else
	{
		while((Hp = &m_lpIMP[Ho>>NA]) < End)
		{
			t	= *Hp;					/* Get IR sample */
			t	*= *Xp;					/* Mult coeff by input sample */
			if(t & 1<<(NHXN-1))			/* Round, if needed */
				t	+= 1<<(NHXN-1);
			t	>>= NHXN;				/* Leave some guard bits, but come back some */
			v	+= t;					/* The filter output */
			Ho	+= m_nDHB;				/* IR step */
			Xp	+= nInc * m_nChannels;	/* Input signal step. NO CHECK ON BOUNDS */
		}
	}
	return v;
}

#endif
