//====================================================================
// File Name : Mixer_Resample.h
// Function  : Audio Mixer Resample interface Declaration file 
// Program   : Li, Shaoen (LSN)
// Date      : Jun 30, 2008
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#ifndef __MIXER_RESAMPLE_H_2008_06_30__
#define __MIXER_RESAMPLE_H_2008_06_30__

#if	1

//====================================================================
//	Audio Mixer Resample Interface class declaration
//====================================================================
class Mixer_Resample
{
public:
	Mixer_Resample()	{}
	virtual	~Mixer_Resample()	{}

	//	Set Sample Method
	virtual	BOOL	SetMethod(UINT nMethod, UINT nChannels) = 0;

	//	Set Sample Ratio
	virtual	BOOL	SetRatio(UINT nOutSampleRate, UINT nInSampleRate) = 0;
	virtual	BOOL	SetRatio(DOUBLE dblRatio) = 0;

	//	Convert
	virtual	INT		ConvertUp(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum) = 0;
	virtual	INT		ConvertDown(SHORT* lpDes, INT nDesFrameNum, SHORT* lpSrc, INT nSrcFrameNum, INT* pUsedFrameNum) = 0;

private:
	//	Copy Construtor
	Mixer_Resample(const Mixer_Resample&);
	//	Operator = overload
	Mixer_Resample&	operator = (const Mixer_Resample&);
};

typedef	Mixer_Resample	MIXER_RESAMPLE;
//====================================================================
//	Audio Mixer Resample Interface class declaration end
//====================================================================

#endif

#endif
