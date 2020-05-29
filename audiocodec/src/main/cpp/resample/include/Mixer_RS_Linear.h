//====================================================================
// File Name : Mixer_RS_Linear.h
// Function  : Audio Mixer Resample Linear Declaration file 
// Program   : Li, Shaoen (LSN)
// Date      : Jun 30, 2008
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#ifndef __MIXER_RS_LINEAR_H_2008_06_30__
#define __MIXER_RS_LINEAR_H_2008_06_30__

#if	MIXER_API_EN

//====================================================================
//	Audio Mixer Resample constant
//====================================================================
//	Conversion constants
#define		NHC			8
#define		NA			7
#define		NP			(NHC+NA)
#define		NPC			(1<<NHC)
#define		AMASK		((1<<NA)-1)
#define		PMASK		((1<<NP)-1)
#define		NH			16
#define		NB			16
#define		NHXN		14
#define		NHG			(NH-NHXN)
#define		NLPSCL		13

#define		MAX_SHORT	(32767)
#define		MIN_SHORT	(-32768)

#ifndef	PI
#define		PI			(3.14159265358979232846)
#endif

#ifndef	PI2
#define		PI2			(6.28318530717958465692)
#endif

#define		D2R			(0.01745329348)			/* (2*pi)/360 */
#define		R2D			(57.29577951)			/* 360/(2*pi) */

#ifndef	MAX
#define		MAX(x,y)	((x)>(y) ?(x):(y))
#endif
#ifndef	MIN
#define		MIN(x,y)	((x)<(y) ?(x):(y))
#endif

#ifndef	ABS
#define		ABS(x)		((x)<0 ? (-(x)) : (x))
#endif

#ifndef	SGN
#define		SGN(x)		((x)<0 ? (-1) : ((x)==0?(0):(1)))
#endif

#define		MAX_CHANNEL	2
//====================================================================
//	Audio Mixer Resample constant End
//====================================================================


//====================================================================
//	Audio Mixer Resample Linear class declaration
//====================================================================
class Mixer_RS_Linear : public Mixer_Resample
{
public:
	Mixer_RS_Linear();
	virtual	~Mixer_RS_Linear();

	//	Set Sample Method
	virtual	BOOL	SetMethod(UINT nMethod, UINT nChannels);

	//	Set Sample Ratio
	virtual	BOOL	SetRatio(UINT nOutSampleRate, UINT nInSampleRate);
	virtual	BOOL	SetRatio(DOUBLE dblRatio);

	// Reset
	virtual	BOOL	Reset();

	//	Convert
	virtual	INT		ConvertUp(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum);
	virtual	INT		ConvertDown(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum);

protected:
	//	Int to Short
	SHORT	IntToShort(INT nVal);

private:
	//	Copy Construtor
	Mixer_RS_Linear(const Mixer_RS_Linear&);
	//	Operator = overload
	Mixer_RS_Linear&	operator = (const Mixer_RS_Linear&);

protected:
	//	Channels
	INT		m_nChannels;
	//	Output sample step
	INT		m_nDTB;
	//	Last Sample
	INT		m_nLastSample[MAX_CHANNEL];
	//	Last Ration
	INT		m_nLastRatio;
	//	Next blk index
	INT		m_nNextInx;
	//	Sample Ratio
	DOUBLE	m_dblRatio;
};

typedef	Mixer_RS_Linear	MIXER_RS_LINEAR;
//====================================================================
//	Audio Mixer Resample Linear class declaration end
//====================================================================

#endif

#endif
