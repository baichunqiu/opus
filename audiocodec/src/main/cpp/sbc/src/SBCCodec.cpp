//#include "StdAfx.h"
#include <memory.h>
#include "SBCCodec.h"
#include "sbc.h"


CSBCCodec::CSBCCodec(
	)
{
	memset(&decoder_sbc,0,sizeof(decoder_sbc));
	memset(&encoder_sbc,0,sizeof(encoder_sbc));
}

CSBCCodec::~CSBCCodec()
{
}

void CSBCCodec::sbc_decoder_initialize()
{
	sbc_init(&decoder_sbc, 0L);
	//decoder_info.codesize = sbc_get_codesize(&decoder_sbc);
	decoder_sbc.endian = SBC_LE;
	//decoder_sbc.endian = SBC_BE;
}
//char *filename, int subbands, int bitpool, int joint,
//int dualchannel, int snr, int blocks

int CSBCCodec::sbc_encoder_initialize(int nSampleRate,int nChannels,int subbands, int bitpool, int joint,int dualchannel, int snr, int blocks)
{
	encoder_info.samplerate = nSampleRate;
//	struct sbc_struct {
//	unsigned long flags;
//
//	uint8_t frequency;
//	uint8_t blocks;
//	uint8_t subbands;
//	uint8_t mode;
//	uint8_t allocation;
//	uint8_t bitpool;
//	uint8_t endian;
//
//	void *priv;
//	void *priv_alloc_base;
//};
	//encoder_sbc.samplerate = nSampleRate;

	encoder_info.channels = nChannels;
	encoder_info.bitpool = bitpool;
	//encoder_info.encode_mode = stereomode;
	encoder_info.blocks = blocks;
	sbc_init(&encoder_sbc, 0L);

	switch ((encoder_info.samplerate)) {
	case 16000:
		encoder_sbc.frequency = SBC_FREQ_16000;
		break;
	case 32000:
		encoder_sbc.frequency = SBC_FREQ_32000;
		break;
	case 44100:
		encoder_sbc.frequency = SBC_FREQ_44100;
		break;
	case 48000:
		encoder_sbc.frequency = SBC_FREQ_48000;
		break;
	}

	encoder_sbc.subbands = subbands == 4? SBC_SB_4 : SBC_SB_8;
	//encoder_sbc.subbands = subbands == 4? 4 : 8;
	if ((encoder_info.channels) == 1) {
		encoder_sbc.mode = SBC_MODE_MONO;
		if (joint || dualchannel) {
			return -1;
		}
	} else if (joint && !dualchannel)
		encoder_sbc.mode = SBC_MODE_JOINT_STEREO;
	else if (!joint && dualchannel)
		encoder_sbc.mode = SBC_MODE_DUAL_CHANNEL;
	else if (!joint && !dualchannel)
		encoder_sbc.mode = SBC_MODE_STEREO;
	else {
		return -2;
	}

	/*encoder_sbc.endian = SBC_LE;*/
	encoder_sbc.endian = SBC_LE;

	encoder_sbc.bitpool = bitpool;
	encoder_sbc.allocation = snr ? SBC_AM_SNR : SBC_AM_LOUDNESS;
	encoder_sbc.blocks = blocks == 4 ? SBC_BLK_4 :
			blocks == 8 ? SBC_BLK_8 :
				blocks == 12 ? SBC_BLK_12 : SBC_BLK_16;
	encoder_info.codesize = sbc_get_codesize(&encoder_sbc);

}
int CSBCCodec::encode(const void *input, int input_len,void *output, int output_len,  int *readed,int* encoded )
{
	//int framelen = 0;
	int nframes = 0;
	int len_res = 0;
	int	encoded_pcm_len = 0;	//已经解码为PCM的长度
	int output_buf_size = 0;	//输出buffer的长度
	int	input_buf_size = 0;		//输入buffer的长度
	unsigned int	frame_encoded_pcm_len = 0;	//每次编码的PCM数据长度
	int	frame_encoded_sbc_len = 0;	//每帧压缩以后，SBC的长度
	unsigned int	encoded_sbc_len = 0;	//已经解码的SBC流的长度

	output_buf_size = output_len;
	input_buf_size = input_len;

	*readed = 0;
	*encoded = 0;

	frame_encoded_pcm_len = encoder_info.codesize;
	//according to the framesize;
	if (input_buf_size % frame_encoded_pcm_len != 0)
		return -1;

	nframes = input_buf_size/frame_encoded_pcm_len;

	encoded_pcm_len = 0;
	encoded_sbc_len = 0;

	while (encoded_pcm_len < input_buf_size) {
		len_res = sbc_encode(&encoder_sbc, (char*)input + encoded_pcm_len , frame_encoded_pcm_len,
			(char*)output + encoded_sbc_len, output_buf_size - encoded_sbc_len,	&frame_encoded_sbc_len);
		if (len_res != frame_encoded_pcm_len || frame_encoded_sbc_len <= 0) {
			return -2;
		}
		encoded_pcm_len += len_res;
		encoded_sbc_len += frame_encoded_sbc_len;
	}
	/*while(input_buf_size>=frame_encoded_pcm_len){
		len_res = sbc_encode(&encoder_sbc, (char*)input + encoded_pcm_len , 512,
			(char*)output + encoded_sbc_len, 123*4,	&frame_encoded_sbc_len);
		encoded_pcm_len += len_res;
		encoded_sbc_len += frame_encoded_sbc_len;
	}*/



	*readed = encoded_pcm_len;
	*encoded = encoded_sbc_len;

	return 0;
}
int CSBCCodec::decode(const void *input, int input_len,void *output, int output_len, int *written,int* decoded )
{
	int framelen = 0;
	int	decoded_pcm_len = 0;	//已经解码为PCM的长度
	int output_buf_size = 0;	//输出buffer的长度
	int	input_buf_size = 0;		//输入buffer的长度
	unsigned int frame_decoded_pcm_len = 0;	//每次解码后PCM数据的长度
	unsigned int	decoded_sbc_len = 0;	//已经解码的SBC流的长度

	output_buf_size = output_len;
	input_buf_size = input_len;

	framelen = sbc_decode(&decoder_sbc, (char*)input, input_buf_size, (char*)output, output_buf_size, &frame_decoded_pcm_len);

	//decoder_info.encode_mode = decoder_sbc.mode == SBC_MODE_MONO ? 1 : 2;
	decoder_info.channels = decoder_sbc.mode == SBC_MODE_MONO ? 1 : 2;
	decoder_info.encode_mode = decoder_sbc.mode == SBC_MODE_MONO ? MONO :decoder_sbc.mode == SBC_MODE_STEREO ? STEREO : JOINT_STEREO;
	switch (decoder_sbc.frequency) {
	case SBC_FREQ_16000:
		decoder_info.samplerate = 16000;
		break;
	case SBC_FREQ_32000:
		decoder_info.samplerate = 32000;
		break;
	case SBC_FREQ_44100:
		decoder_info.samplerate = 44100;
		break;
	case SBC_FREQ_48000:
		decoder_info.samplerate = 48000;
		break;
	default:
		decoder_info.samplerate = 0;
		break;
	}

	decoded_sbc_len += framelen;
	decoded_pcm_len += frame_decoded_pcm_len;

	while (framelen > 0) {
		if (decoded_pcm_len > output_buf_size)
			break;

		/* push the pointer in the file forward to the next bit to be
		 * decoded tell the decoder to decode up to the remaining
		 * length of the file (!) */
		
		framelen = sbc_decode(&decoder_sbc, (char*)input + decoded_sbc_len, input_buf_size - decoded_sbc_len,
					(char*)output + decoded_pcm_len, sizeof(output) - decoded_pcm_len, &frame_decoded_pcm_len);
		if (framelen < 0)
			break;
		decoded_sbc_len += framelen;
		decoded_pcm_len += frame_decoded_pcm_len;
	}


	//while (framelen > 0) {
	//	if (decoded_pcm_len > output_buf_size)
	//		break;

	//	/* push the pointer in the file forward to the next bit to be
	//	 * decoded tell the decoder to decode up to the remaining
	//	 * length of the file (!) */
	//	
	//	framelen = sbc_decode(&decoder_sbc, (char*)input + decoded_sbc_len,123,
	//				(char*)output + decoded_pcm_len, 512, &frame_decoded_pcm_len);
	//	if (framelen < 0)
	//		break;
	//	decoded_sbc_len += framelen;
	//	decoded_pcm_len += frame_decoded_pcm_len;
	//}


	*written = decoded_pcm_len;
	*decoded = decoded_sbc_len;
	return frame_decoded_pcm_len;
}

int CSBCCodec::be_to_le(void *input,int nLen)
{
	//the endian swap
	if(nLen % 2 != 0)
		return -1;
	unsigned short *pshort = (unsigned short *)input;
	unsigned short nData = 0;
	for (int i = 0;i < nLen/2;i++)
	{
		nData = *(pshort + i);
		nData = bswap_16(nData);
		*(pshort + i) = nData;
	}
	return 0;
}

sbc_codec_info CSBCCodec::get_encoder_info()
{
	return encoder_info;
}