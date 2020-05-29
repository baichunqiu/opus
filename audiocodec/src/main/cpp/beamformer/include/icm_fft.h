#ifndef ICM_FFT_H
#define ICM_FFT_H

#include "complex.h"

typedef struct fft_t fft_t;
typedef struct cdft_t cdft_t;
typedef struct kft_t kft_t;

#ifdef  __cplusplus
extern "C" {
#endif

fft_t* icm_fft_init(OsInt32 M);
void icm_fft_uninit(fft_t *fft);
// dc 0 re im re im nyquist 0
void icm_fft(fft_t *fft,OsFloat *x,complex *xf);
void icm_ifft(fft_t *fft,complex *xf,OsFloat *x);
// dc nyquist re im re im
void icm_fft2(fft_t *fft,OsFloat *x);
void icm_ifft2(fft_t *fft,OsFloat *x);

// from speex (fft of arbitrary length)
kft_t* kiss_fft_init(OsInt32 M);
void kiss_fft_uninit(kft_t *fft);

void kiss_fft(kft_t *fft,OsFloat *x,complex *fx);
void kiss_ifft(kft_t *fft,complex *fx,OsFloat *x);
void kiss_fft2(kft_t *fft,OsFloat *x);
void kiss_ifft2(kft_t *fft,OsFloat *x);

cdft_t* icm_cdft_init(OsInt32 M);
void icm_cdft_uninit(cdft_t *fft);
void icm_cdft_fft(cdft_t *fft,complex x[]);
void icm_cdft_ifft(cdft_t *fft,complex x[]);

#ifdef  __cplusplus
}
#endif

#endif