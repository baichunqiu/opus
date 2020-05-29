#ifndef BEAMFORMER_H
#define BEAMFORMER_H

#include "Platform.h"

#define kNumMics    2

typedef struct Beamformer Beamformer;

#ifdef  __cplusplus
extern "C" {
#endif

/**
 * @brief           打开波束形成器
 * @inGeometry      阵列的形状（每个阵元的3维位置信息）
 * @inFs            采样率
 * @inFrameSize     帧长采样数
 * @return          返回波束形成器
 */
Beamformer* IcmBeamformerOpen(OsFloat (*inGeometry)[3],OsInt32 inFs,OsInt32 inFrameSize);

/**
 * @brief           设置波束形成器参数
 * @inBeamformer    波束形成器实例
 * @inTargetAngle   期望信号的方向，单位：弧度
 * @inInterfAngle   干扰信号的方向，单位：弧度
 */
void IcmBeamformerSetParamer(Beamformer *inBeamformer,OsFloat inTargetAngle,OsFloat inInterfAngle);

/**
 * @brief           关闭波束形成器
 * @inBeamformer    波束形成器实例
 */
void IcmBeamformerClose(Beamformer *inBeamformer);

/**
 * @brief           执行波束形成
 * @inBeamformer    波束形成器实例
 * @inPcm           阵列录制的PCM数据
 * @inLen           输入波束形成器的PCM数据采样数
 * @ioPcm           波束形成器输出
 * @ioLen           波束形成器输出PCM数据采样数
 */
void IcmBeamformerProcess(Beamformer *inBeamformer,OsInt16 *inPcm,OsInt32 inLen,OsInt16 *ioPcm,OsInt32 *ioLen);

#ifdef  __cplusplus
}
#endif

#endif