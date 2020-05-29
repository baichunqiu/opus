//====================================================================
// File Name : Mixer_Def.h
// Function  : Mixer Definition file 
// Program   : Li, Shaoen (LSN)
// Date      : Jan 30, 2007
// Version   : 0.1
// History
//		Version		Author		Memo	
//		0.1			Lishaoen	Programming start
//====================================================================

#ifndef __MIXER_DEF_H_2007_01_30__
#define __MIXER_DEF_H_2007_01_30__


#include "Def_DataType.h"

#define MIXER_API_EN   1
//====================================================================
//	Mixer Constant definition
//====================================================================
//	EQ Mode definition
#define		EQ_MODE_NORMAL		1
#define		EQ_MODE_CLASSIC		2
#define		EQ_MODE_JAZZ		3
#define		EQ_MODE_POP			4
#define		EQ_MODE_ROCK		5
//====================================================================
//	Mixer Constant definition End
//====================================================================

//====================================================================
//	Mixer Volume definition
//====================================================================
#define		VOLUME_MIN			0
#define		VOLUME_MAX			0x7FFF
#define		VOLUME_DEF			0x100
//====================================================================
//	Mixer Volume definition End
//====================================================================

//====================================================================
//	Mixer Buffer definition
//====================================================================
#define		DEF_REPLAY_BUFLEN	16*1024
#define		DEF_MIXER_BUFLEN	16*1024
#define		DEF_RECORD_BUFLEN	16*1024
//====================================================================
//	Mixer Buffer definition End
//====================================================================

//====================================================================
//	Mixer Channel definition
//====================================================================
#define		CHANNEL_MONO		0x00000001
#define		CHANNEL_STEREO		0x00000002
//====================================================================
//	Mixer Channel definition End
//====================================================================

//====================================================================
//	Mixer IOCtl definition
//====================================================================
#define		MIXER_CMD_MASK				0x0000000F
#define		MIXER_CMD_GETVOLUME			0x00000001
#define		MIXER_CMD_SETVOLUME			0x00000002
#define		MIXER_CMD_GETMUTE			0x00000003
#define		MIXER_CMD_SETMUTE			0x00000004
#define		MIXER_CMD_GETCHANNEL		0x00000005
#define		MIXER_CMD_SETCHANNEL		0x00000006
#define		MIXER_CMD_GETEQ				0x00000007
#define		MIXER_CMD_SETEQ				0x00000008
#define		MIXER_CMD_STOPALL			0x00000009
//====================================================================
//	Mixer IOCtl definition End
//====================================================================

//====================================================================
//	Mixer type definition
//====================================================================
#define		DEF_SAMPLE_LEN		16
#define		DEF_SAMPLE_RATE		48000

//	Audio FMT Info
typedef	struct _tagAUDIOFMTINFO
{
	UINT	m_nSampleLen;
	UINT	m_nSampleRate;
	UINT	m_nChannel;
	UINT	m_nEQMode;
}AUDIOFMTINFO, *LPAUDIOFMTINFO;

//	Audio Play Status
#define		PS_STOP				0x00000000
#define		PS_PLAY				0x00000001
#define		PS_SUSPEND			0x00000002

//	Audio Ctrl Info
typedef	struct _tagAUDIOCTRLINFO
{
	BOOL	m_bMute;
	BOOL	m_nSuspend;
	UINT	m_nVolume;
	UINT	m_nStatus;
	INT		m_nCurPlayPos;
	INT		m_nCurFillLen;
}AUDIOCTRLINFO, *LPAUDIOCTRLINFO;

//	Audio Capacity Info
typedef	struct _tagAUDIOCAPACITY
{
	LPVOID	m_lpBuf;
	UINT	m_nLen;
	INT		m_nCurPlayPos;
	INT		m_nCurFillLen;
}AUDIOCAPACITY, *LPAUDIOCAPACITY;
//====================================================================
//	Mixer type definition End
//====================================================================

//====================================================================
//	Mixer Error Code definition
//====================================================================
typedef		INT		AS_STATUS;

#define		AS_OK					0
#define		AS_ERR					-1
#define		AS_ERR_INVALIDHANDLE	1
#define		AS_ERR_INVALIDPARAM		2
#define		AS_ERR_INVALIDBUFFER	3
#define		AS_ERR_MEMALLOCATEFAIL	4
#define		AS_ERR_CHANNELSUSPEND	5
#define		AS_ERR_BUFFERISFULL		6
#define		AS_ERR_ISPLAYING		7
#define		AS_ERR_ISNOTPLAYING		8

//====================================================================
//	Mixer Error Code definition End
//====================================================================

#ifdef	__cplusplus
//====================================================================
//	Mixer Data Type definition
//====================================================================
#include "ArrayTemplate.h"
typedef		CArray1D<BYTE>		BYTE_BUFFER;
typedef		CArray1D<SHORT>		PCM_BUFFER;
//====================================================================
//	Mixer Data Type definition End
//====================================================================
#endif

#endif
