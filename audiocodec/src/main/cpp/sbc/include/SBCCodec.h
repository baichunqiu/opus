#pragma once
#include "sbc.h"

typedef enum 
{
	MONO,
	JOINT_STEREO,
	STEREO
}SBC_ENCODE_MODE;

typedef enum 
{
	SNR,
	LOUDNESS
}SBC_ALLOC_MODE;

typedef struct
{
	int samplerate;
	int subbands;
	int	bitpool;
	int channels;
	int	blocks;
	int	codesize;
	SBC_ENCODE_MODE encode_mode;
	SBC_ALLOC_MODE	aloc_mode;
}sbc_codec_info;


#define		CODEC_BUF_SIZE	8192

class CSBCCodec
{
public:
	
public:
	CSBCCodec();
	~CSBCCodec();
public:
	void sbc_decoder_initialize();
	int sbc_encoder_initialize(int nSampleRate,int nChannels,int subbands, int bitpool, int joint,int dualchannel, int snr, int blocks);
	sbc_codec_info get_decoder_info();
	sbc_codec_info get_encoder_info();

	//return the frame size
	//input is the input stream buffer
	//output is the internal buffer,outputlen is the 
	int decode(const void *input, int input_len,void *output, int output_len,  int *written,int* decoded );
	int encode(const void *input, int input_len,void *output, int output_len,  int *written,int* decoded );

	int	be_to_le(void *input,int nLen);

private:
	sbc_t	encoder_sbc;
	sbc_t	decoder_sbc;

	int		decoder_frame_len;
	int		decoded_frame_len;

	sbc_codec_info	decoder_info;
	sbc_codec_info	encoder_info;

private:
	char	decoder_buf[CODEC_BUF_SIZE];
	char	encoder_buf[CODEC_BUF_SIZE];
};

