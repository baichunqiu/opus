//====================================================================
// File Name : Mixer_RS_KW.h
// Function  : Audio Mixer Resample Kaiser-Windowed Declaration file 
// Program   : Li, Shaoen (LSN)
// Date      : Jun 30, 2008
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#ifndef __MIXER_RS_KW_H_2008_06_30__
#define __MIXER_RS_KW_H_2008_06_30__

#if	1

//====================================================================
//	Audio Mixer Resample Kaiser-Windowed class declaration
//====================================================================
class Mixer_RS_KW : public Mixer_RS_Linear
{
public:
	Mixer_RS_KW();
	virtual	~Mixer_RS_KW();

	//	Set Sample Method
	virtual	BOOL	SetMethod(UINT nMethod, UINT nChannels);

	//	Set Sample Ratio
	virtual	BOOL	SetRatio(UINT nOutSampleRate, UINT nInSampleRate);
	virtual	BOOL	SetRatio(DOUBLE dblRatio);

	//	Set Buffer Length
			BOOL	SetPCMBufLen(UINT nBufLen);

	//	Reset
			BOOL	Reset();

	//	Convert
	virtual	INT		ConvertUp(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum);
	virtual	INT		ConvertDown(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum);

protected:
	//	Int to Short
	SHORT	IntToShort(INT nVal);
	//	Filter UP
	INT		FilterUp(SHORT* lpSrc, INT nPhase, INT nInc);
	//	Filter Down
	INT		FilterDown(SHORT* lpSrc, INT nPhase, INT nInc);

private:
	//	Copy Construtor
	Mixer_RS_KW(const Mixer_RS_KW&);
	//	Operator = overload
	Mixer_RS_KW&	operator = (const Mixer_RS_KW&);

protected:
	//	The Impulse response, coefficients
	SHORT*	m_lpIMP;
	//	The Impulse response deltas, ImpD[n] = Imp[n+1]-Imp[n]
	SHORT*	m_lpIMPD;
	//	The Filter table length
	INT		m_nWing;
	//	The Normalize for unity filter gain
	INT		m_nSCL;
	//	The Mult
	INT		m_nMult;
	//	The DHB
	INT		m_nDHB;
	//	The X Offset
	INT		m_nXOff;
	//	The Interpolation
	BOOL	m_bInterp;
	//	The Local Buffer
	PCM_BUFFER	m_aryPCM;
};

typedef	Mixer_RS_KW	MIXER_RS_KW;
//====================================================================
//	Audio Mixer Resample Kaiser-Windowed class declaration end
//====================================================================

#endif

#endif
