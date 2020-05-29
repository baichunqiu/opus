#include "../include/icm_fft.h"

void rdft(OsInt32,OsInt32,OsFloat *,OsInt32 *,OsFloat *);
void cdft(OsInt32,OsInt32,OsFloat *,OsInt32 *,OsFloat *);

struct fft_t
{
    OsInt32 M;
    OsInt32 freq_bins;

    OsFloat *x;

    OsInt32 *ip;    // work area for bit reversal L = M
    OsFloat *w;     // cos/sin table, are initialized if ip[0] == 0. L = M/2
};

fft_t* icm_fft_init(OsInt32 M)
{
    OsInt32 bits = ceil(log(1.0*M)/log(2.0));
    fft_t *fft = 0;

    if(M != pow(2.0,bits)) return 0;

    fft = (fft_t*)icm_alloc(1,sizeof(fft_t));
    assert(0 != fft);
    fft->M = M;
    fft->freq_bins = M/2+1;

    fft->x = (OsFloat*)icm_alloc(M,sizeof(OsFloat));
    assert(0 != fft->x);
    fft->ip = (OsInt32*)icm_alloc(M,sizeof(OsInt32));
    assert(0 != fft->ip);
    fft->w  = (OsFloat*)icm_alloc(M/2,sizeof(OsFloat));
    assert(0 != fft->w);

    fft->ip[0] = 0;
    rdft(M,1,fft->x,fft->ip,fft->w);
    return fft;
}

void icm_fft_uninit(fft_t *fft)
{
    if(fft->x) icm_free(fft->x);
    icm_free(fft->ip);
    icm_free(fft->w);
    icm_free(fft);
}

void icm_fft(fft_t *fft,OsFloat *x,complex *fx)
{
    OsInt32 i;

    memcpy(fft->x,x,fft->M*sizeof(OsFloat));
    rdft(fft->M,1,fft->x,fft->ip,fft->w);

    fx[0] = init_c(fft->x[0],0);
    fx[fft->M/2] = init_c(fft->x[1],0);
    for(i = 1; i < fft->M/2; ++i)
    {
        fx[i] = init_c(fft->x[2*i],-fft->x[2*i+1]);
    }
}

void icm_ifft(fft_t *fft,complex *fx,OsFloat *x)
{
    OsInt32 N = fft->M/2;
    OsFloat scale = 2.0/fft->M;
    OsInt32 i;

    x[0] = fx[0].re;
    x[1] = fx[N].re;
    for(i = 1; i < N; ++i)
    {
        x[2*i]  = fx[i].re;
        x[2*i+1]= -fx[i].im;
    }
    rdft(fft->M,-1,x,fft->ip,fft->w);

    for(i = 0; i < fft->M; i++)
    {
        x[i] *= scale;
    }
}

void icm_fft2(fft_t *fft,OsFloat *x)
{
    OsInt32 i;

    rdft(fft->M,1,x,fft->ip,fft->w);
    for(i = 3; i < fft->M; i += 2)
    {
        x[i] = -x[i];
    }
}

void icm_ifft2(fft_t *fft,OsFloat *x)
{
    OsInt32 N = fft->M/2;
    OsInt32 i;

    for(i = 3; i < fft->M; i += 2)
    {
        x[i] = -x[i];
    }
    rdft(fft->M,-1,x,fft->ip,fft->w);
    for(i = 0; i < fft->M; i++)
    {
        x[i] /= N;
    }
}

/*
Fast Fourier/Cosine/Sine Transform
    dimension   :one
    data length :power of 2
    decimation  :frequency
    radix       :4, 2
    data        :inplace
    table       :use
functions
    cdft: Complex Discrete Fourier Transform
    rdft: Real Discrete Fourier Transform
    ddct: Discrete Cosine Transform
    ddst: Discrete Sine Transform
    dfct: Cosine Transform of RDFT (Real Symmetric DFT)
    dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
function prototypes
    void cdft(OsInt32, OsInt32, OsFloat *, OsInt32 *, OsFloat *);
    void rdft(OsInt32, OsInt32, OsFloat *, OsInt32 *, OsFloat *);
    void ddct(OsInt32, OsInt32, OsFloat *, OsInt32 *, OsFloat *);
    void ddst(OsInt32, OsInt32, OsFloat *, OsInt32 *, OsFloat *);
    void dfct(OsInt32, OsFloat *, OsFloat *, OsInt32 *, OsFloat *);
    void dfst(OsInt32, OsFloat *, OsFloat *, OsInt32 *, OsFloat *);


-------- Complex DFT (Discrete Fourier Transform) --------
    [definition]
        <case1>
            X[k] = sum_j=0^n-1 x[j]*exp(2*pi*i*j*k/n), 0<=k<n
        <case2>
            X[k] = sum_j=0^n-1 x[j]*exp(-2*pi*i*j*k/n), 0<=k<n
        (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
    [usage]
        <case1>
            ip[0] = 0; // first time only
            cdft(2*n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            cdft(2*n, -1, a, ip, w);
    [parameters]
        2*n            :data length (OsInt32)
                        n >= 1, n = power of 2
        a[0...2*n-1]   :input/output data (OsFloat *)
                        input data
                            a[2*j] = Re(x[j]), 
                            a[2*j+1] = Im(x[j]), 0<=j<n
                        output data
                            a[2*k] = Re(X[k]), 
                            a[2*k+1] = Im(X[k]), 0<=k<n
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            cdft(2*n, -1, a, ip, w);
        is 
            cdft(2*n, 1, a, ip, w);
            for (j = 0; j <= 2 * n - 1; j++) {
                a[j] *= 1.0 / n;
            }
        .


-------- Real DFT / Inverse of Real DFT --------
    [definition]
        <case1> RDFT
            R[k] = sum_j=0^n-1 a[j]*cos(2*pi*j*k/n), 0<=k<=n/2
            I[k] = sum_j=0^n-1 a[j]*sin(2*pi*j*k/n), 0<k<n/2
        <case2> IRDFT (excluding scale)
            a[k] = (R[0] + R[n/2]*cos(pi*k))/2 + 
                   sum_j=1^n/2-1 R[j]*cos(2*pi*j*k/n) + 
                   sum_j=1^n/2-1 I[j]*sin(2*pi*j*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            rdft(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            rdft(n, -1, a, ip, w);
    [parameters]
        n              :data length (OsInt32)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (OsFloat *)
                        <case1>
                            output data
                                a[2*k] = R[k], 0<=k<n/2
                                a[2*k+1] = I[k], 0<k<n/2
                                a[1] = R[n/2]
                        <case2>
                            input data
                                a[2*j] = R[j], 0<=j<n/2
                                a[2*j+1] = I[j], 0<j<n/2
                                a[1] = R[n/2]
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n/2-1]   :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            rdft(n, 1, a, ip, w);
        is 
            rdft(n, -1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
    [definition]
        <case1> IDCT (excluding scale)
            C[k] = sum_j=0^n-1 a[j]*cos(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DCT
            C[k] = sum_j=0^n-1 a[j]*cos(pi*(j+1/2)*k/n), 0<=k<n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddct(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddct(n, -1, a, ip, w);
    [parameters]
        n              :data length (OsInt32)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (OsFloat *)
                        output data
                            a[k] = C[k], 0<=k<n
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddct(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddct(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- DST (Discrete Sine Transform) / Inverse of DST --------
    [definition]
        <case1> IDST (excluding scale)
            S[k] = sum_j=1^n A[j]*sin(pi*j*(k+1/2)/n), 0<=k<n
        <case2> DST
            S[k] = sum_j=0^n-1 a[j]*sin(pi*(j+1/2)*k/n), 0<k<=n
    [usage]
        <case1>
            ip[0] = 0; // first time only
            ddst(n, 1, a, ip, w);
        <case2>
            ip[0] = 0; // first time only
            ddst(n, -1, a, ip, w);
    [parameters]
        n              :data length (OsInt32)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (OsFloat *)
                        <case1>
                            input data
                                a[j] = A[j], 0<j<n
                                a[0] = A[n]
                            output data
                                a[k] = S[k], 0<=k<n
                        <case2>
                            output data
                                a[k] = S[k], 0<k<n
                                a[0] = S[n]
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n/2)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n/2+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/4-1] :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            ddst(n, -1, a, ip, w);
        is 
            a[0] *= 0.5;
            ddst(n, 1, a, ip, w);
            for (j = 0; j <= n - 1; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Cosine Transform of RDFT (Real Symmetric DFT) --------
    [definition]
        C[k] = sum_j=0^n a[j]*cos(pi*j*k/n), 0<=k<=n
    [usage]
        ip[0] = 0; // first time only
        dfct(n, a, t, ip, w);
    [parameters]
        n              :data length - 1 (OsInt32)
                        n >= 2, n = power of 2
        a[0...n]       :input/output data (OsFloat *)
                        output data
                            a[k] = C[k], 0<=k<=n
        t[0...n/2]     :work area (OsFloat *)
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
        is 
            a[0] *= 0.5;
            a[n] *= 0.5;
            dfct(n, a, t, ip, w);
            for (j = 0; j <= n; j++) {
                a[j] *= 2.0 / n;
            }
        .


-------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
    [definition]
        S[k] = sum_j=1^n-1 a[j]*sin(pi*j*k/n), 0<k<n
    [usage]
        ip[0] = 0; // first time only
        dfst(n, a, t, ip, w);
    [parameters]
        n              :data length + 1 (OsInt32)
                        n >= 2, n = power of 2
        a[0...n-1]     :input/output data (OsFloat *)
                        output data
                            a[k] = S[k], 0<k<n
                        (a[0] is used for work area)
        t[0...n/2-1]   :work area (OsFloat *)
        ip[0...*]      :work area for bit reversal (OsInt32 *)
                        length of ip >= 2+sqrt(n/4)
                        strictly, 
                        length of ip >= 2+(1<<(OsInt32)(log(n/4+0.5)/log(2))/2).
                        ip[0],ip[1] are pointers of the cos/sin table.
        w[0...n*5/8-1] :cos/sin table (OsFloat *)
                        w[],ip[] are initialized if ip[0] == 0.
    [remark]
        Inverse of 
            dfst(n, a, t, ip, w);
        is 
            dfst(n, a, t, ip, w);
            for (j = 1; j <= n - 1; j++)
            {
                a[j] *= 2.0 / n;
            }
Appendix :
    The cos/sin table is recalculated when the larger table required.
    w[] and ip[] are compatible with all routines.
*/


void cdft(OsInt32 n, OsInt32 isgn, OsFloat *a, OsInt32 *ip, OsFloat *w)
{
    void makewt(OsInt32 nw, OsInt32 *ip, OsFloat *w);
    void bitrv2(OsInt32 n, OsInt32 *ip, OsFloat *a);
    void bitrv2conj(OsInt32 n, OsInt32 *ip, OsFloat *a);
    void cftfsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void cftbsub(OsInt32 n, OsFloat *a, OsFloat *w);
    
    if (n > (ip[0] << 2))
    {
        makewt(n >> 2, ip, w);
    }
    if (n > 4)
    {
        if (isgn >= 0)
        {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
        }
        else
        {
            bitrv2conj(n, ip + 2, a);
            cftbsub(n, a, w);
        }
    }
    else if (n == 4)
    {
        cftfsub(n, a, w);
    }
}


void rdft(OsInt32 n,OsInt32 isgn,OsFloat *a,OsInt32 *ip,OsFloat *w)
{
    void makewt(OsInt32 nw,OsInt32 *ip,OsFloat *w);
    void makect(OsInt32 nc,OsInt32 *ip,OsFloat *c);
    void bitrv2(OsInt32 n,OsInt32 *ip,OsFloat *a);
    void cftfsub(OsInt32 n,OsFloat *a,OsFloat *w);
    void cftbsub(OsInt32 n,OsFloat *a,OsFloat *w);
    void rftfsub(OsInt32 n,OsFloat *a,OsInt32 nc,OsFloat *c);
    void rftbsub(OsInt32 n,OsFloat *a,OsInt32 nc,OsFloat *c);
    OsInt32 nw,nc;
    OsFloat xi;
    
    nw = ip[0];
    if (n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw,ip,w);
    }
    nc = ip[1];
    if (n > (nc << 2))
    {
        nc = n >> 2;
        makect(nc,ip,w+nw);
    }
    if (isgn >= 0)
    {
        if (n > 4)
        {
            bitrv2(n,ip+2,a);
            cftfsub(n,a,w);
            rftfsub(n,a,nc,w+nw);
        }
        else if (n == 4)
        {
            cftfsub(n,a,w);
        }
        xi = a[0]-a[1];
        a[0] += a[1];
        a[1] = xi;
    }
    else
    {
        a[1] = 0.5f * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4)
        {
            rftbsub(n,a,nc,w+nw);
            bitrv2(n,ip+2,a);
            cftbsub(n,a,w);
        }
        else if (n == 4)
        {
            cftfsub(n,a,w);
        }
    }
}

void ddct(OsInt32 n, OsInt32 isgn, OsFloat *a, OsInt32 *ip, OsFloat *w)
{
    void makewt(OsInt32 nw, OsInt32 *ip, OsFloat *w);
    void makect(OsInt32 nc, OsInt32 *ip, OsFloat *c);
    void bitrv2(OsInt32 n, OsInt32 *ip, OsFloat *a);
    void cftfsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void cftbsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void rftfsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    void rftbsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    void dctsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    OsInt32 j, nw, nc;
    OsFloat xr;
    
    nw = ip[0];
    if (n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc)
    {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0)
    {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2)
        {
            a[j + 1] = a[j] - a[j - 1];
            a[j] += a[j - 1];
        }
        a[1] = a[0] - xr;
        a[0] += xr;
        if (n > 4)
        {
            rftbsub(n, a, nc, w + nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        }
        else if (n == 4)
        {
            cftfsub(n, a, w);
        }
    }
    dctsub(n,a,nc,w + nw);
    if (isgn >= 0)
    {
        if (n > 4)
        {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        }
        else if (n == 4)
        {
            cftfsub(n, a, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2)
        {
            a[j - 1] = a[j] - a[j + 1];
            a[j] += a[j + 1];
        }
        a[n - 1] = xr;
    }
}


void ddst(OsInt32 n, OsInt32 isgn, OsFloat *a, OsInt32 *ip, OsFloat *w)
{
    void makewt(OsInt32 nw, OsInt32 *ip, OsFloat *w);
    void makect(OsInt32 nc, OsInt32 *ip, OsFloat *c);
    void bitrv2(OsInt32 n, OsInt32 *ip, OsFloat *a);
    void cftfsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void cftbsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void rftfsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    void rftbsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    void dstsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    OsInt32 j, nw, nc;
    OsFloat xr;
    
    nw = ip[0];
    if (n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc)
    {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0)
    {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2)
        {
            a[j + 1] = -a[j] - a[j - 1];
            a[j] -= a[j - 1];
        }
        a[1] = a[0] + xr;
        a[0] -= xr;
        if (n > 4)
        {
            rftbsub(n,a,nc,w+nw);
            bitrv2(n,ip+2,a);
            cftbsub(n,a,w);
        }
        else if (n == 4)
        {
            cftfsub(n,a,w);
        }
    }
    dstsub(n,a,nc,w+nw);
    if (isgn >= 0)
    {
        if (n > 4)
        {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        }
        else if (n == 4)
        {
            cftfsub(n, a, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2)
        {
            a[j - 1] = -a[j] - a[j + 1];
            a[j] -= a[j + 1];
        }
        a[n - 1] = -xr;
    }
}

void dfct(OsInt32 n, OsFloat *a, OsFloat *t, OsInt32 *ip, OsFloat *w)
{
    void makewt(OsInt32 nw, OsInt32 *ip, OsFloat *w);
    void makect(OsInt32 nc, OsInt32 *ip, OsFloat *c);
    void bitrv2(OsInt32 n, OsInt32 *ip, OsFloat *a);
    void cftfsub(OsInt32 n, OsFloat *a, OsFloat *w);
    void rftfsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    void dctsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c);
    OsInt32 j, k, l, m, mh, nw, nc;
    OsFloat xr, xi, yr, yi;
    
    nw = ip[0];
    if (n > (nw << 3))
    {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1))
    {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    m = n >> 1;
    yi = a[m];
    xi = a[0] + a[n];
    a[0] -= a[n];
    t[0] = xi - yi;
    t[m] = xi + yi;
    if (n > 2)
    {
        mh = m >> 1;
        for (j = 1; j < mh; j++)
        {
            k = m - j;
            xr = a[j] - a[n - j];
            xi = a[j] + a[n - j];
            yr = a[k] - a[n - k];
            yi = a[k] + a[n - k];
            a[j] = xr;
            a[k] = yr;
            t[j] = xi - yi;
            t[k] = xi + yi;
        }
        t[mh] = a[mh] + a[n - mh];
        a[mh] -= a[n - mh];
        dctsub(m, a, nc, w + nw);
        if (m > 4)
        {
            bitrv2(m, ip + 2, a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w + nw);
        }
        else if (m == 4)
        {
            cftfsub(m, a, w);
        }
        a[n - 1] = a[0] - a[1];
        a[1] = a[0] + a[1];
        for (j = m - 2; j >= 2; j -= 2)
        {
            a[2 * j + 1] = a[j] + a[j + 1];
            a[2 * j - 1] = a[j] - a[j + 1];
        }
        l = 2;
        m = mh;
        while (m >= 2)
        {
            dctsub(m, t, nc, w + nw);
            if (m > 4)
            {
                bitrv2(m, ip + 2, t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w + nw);
            }
            else if (m == 4)
            {
                cftfsub(m, t, w);
            }
            a[n - l] = t[0] - t[1];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2)
            {
                k += l << 2;
                a[k - l] = t[j] - t[j + 1];
                a[k + l] = t[j] + t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 0; j < mh; j++)
            {
                k = m - j;
                t[j] = t[m + k] - t[m + j];
                t[k] = t[m + k] + t[m + j];
            }
            t[mh] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
        a[n] = t[2] - t[1];
        a[0] = t[2] + t[1];
    }
    else
    {
        a[1] = a[0];
        a[2] = t[0];
        a[0] = t[1];
    }
}

void dfst(OsInt32 n,OsFloat *a,OsFloat *t,OsInt32 *ip,OsFloat *w)
{
    void makewt(OsInt32 nw,OsInt32 *ip,OsFloat *w);
    void makect(OsInt32 nc,OsInt32 *ip,OsFloat *c);
    void bitrv2(OsInt32 n,OsInt32 *ip,OsFloat *a);
    void cftfsub(OsInt32 n,OsFloat *a,OsFloat *w);
    void rftfsub(OsInt32 n,OsFloat *a,OsInt32 nc,OsFloat *c);
    void dstsub(OsInt32 n,OsFloat *a,OsInt32 nc,OsFloat *c);
    OsInt32 j,k,l,m,mh,nw,nc;
    OsFloat xr,xi,yr,yi;
    
    nw = ip[0];
    if (n > (nw << 3))
    {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1))
    {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    if (n > 2)
    {
        m = n >> 1;
        mh = m >> 1;
        for (j = 1; j < mh; j++)
        {
            k = m - j;
            xr = a[j] + a[n - j];
            xi = a[j] - a[n - j];
            yr = a[k] + a[n - k];
            yi = a[k] - a[n - k];
            a[j] = xr;
            a[k] = yr;
            t[j] = xi + yi;
            t[k] = xi - yi;
        }
        t[0] = a[mh] - a[n - mh];
        a[mh] += a[n - mh];
        a[0] = a[m];
        dstsub(m, a, nc, w + nw);
        if (m > 4)
        {
            bitrv2(m, ip + 2, a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w + nw);
        }
        else if (m == 4)
        {
            cftfsub(m, a, w);
        }
        a[n - 1] = a[1] - a[0];
        a[1] = a[0] + a[1];
        for (j = m - 2; j >= 2; j -= 2)
        {
            a[2 * j + 1] = a[j] - a[j + 1];
            a[2 * j - 1] = -a[j] - a[j + 1];
        }
        l = 2;
        m = mh;
        while (m >= 2)
        {
            dstsub(m, t, nc, w + nw);
            if (m > 4)
            {
                bitrv2(m, ip + 2, t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w + nw);
            }
            else if (m == 4)
            {
                cftfsub(m, t, w);
            }
            a[n - l] = t[1] - t[0];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2)
            {
                k += l << 2;
                a[k - l] = -t[j] - t[j + 1];
                a[k + l] = t[j] - t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 1; j < mh; j++)
            {
                k = m - j;
                t[j] = t[m + k] + t[m + j];
                t[k] = t[m + k] - t[m + j];
            }
            t[0] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
    }
    a[0] = 0;
}

/* -------- initializing routines -------- */

#include <math.h>

void makewt(OsInt32 nw,OsInt32 *ip,OsFloat *w)
{
    void bitrv2(OsInt32 n,OsInt32 *ip,OsFloat *a);
    OsInt32 j, nwh;
    OsFloat delta, x, y;
    
    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2)
    {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        if (nwh > 2)
        {
            for (j = 2; j < nwh; j += 2)
            {
                x = cos(delta * j);
                y = sin(delta * j);
                w[j] = x;
                w[j + 1] = y;
                w[nw - j] = y;
                w[nw - j + 1] = x;
            }
            bitrv2(nw, ip + 2, w);
        }
    }
}

void makect(OsInt32 nc,OsInt32 *ip,OsFloat *c)
{
    OsInt32 j, nch;
    OsFloat delta;
    
    ip[1] = nc;
    if (nc > 1)
    {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = cos(delta * nch);
        c[nch] = 0.5 * c[0];
        for (j = 1; j < nch; j++)
        {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}

/* -------- child routines -------- */

void bitrv2(OsInt32 n,OsInt32 *ip,OsFloat *a)
{
    OsInt32 j, j1, k, k1, l, m, m2;
    OsFloat xr, xi, yr, yi;
    
    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l)
    {
        l >>= 1;
        for (j = 0; j < m; j++)
        {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    m2 = 2 * m;
    if ((m << 3) == l)
    {
        for (k = 0; k < m; k++)
        {
            for (j = 0; j < k; j++)
            {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 -= m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            j1 = 2 * k + m2 + ip[k];
            k1 = j1 + m2;
            xr = a[j1];
            xi = a[j1 + 1];
            yr = a[k1];
            yi = a[k1 + 1];
            a[j1] = yr;
            a[j1 + 1] = yi;
            a[k1] = xr;
            a[k1 + 1] = xi;
        }
    }
    else
    {
        for (k = 1; k < m; k++)
        {
            for (j = 0; j < k; j++)
            {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}

void bitrv2conj(OsInt32 n, OsInt32 *ip,OsFloat *a)
{
    OsInt32 j, j1, k, k1, l, m, m2;
    OsFloat xr, xi, yr, yi;
    
    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l)
    {
        l >>= 1;
        for (j = 0; j < m; j++)
        {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    m2 = 2 * m;
    if ((m << 3) == l)
    {
        for (k = 0; k < m; k++)
        {
            for (j = 0; j < k; j++)
            {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 -= m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            j1 = k1 + m2;
            k1 = j1 + m2;
            xr = a[j1];
            xi = -a[j1 + 1];
            yr = a[k1];
            yi = -a[k1 + 1];
            a[j1] = yr;
            a[j1 + 1] = yi;
            a[k1] = xr;
            a[k1 + 1] = xi;
            k1 += m2;
            a[k1 + 1] = -a[k1 + 1];
        }
    }
    else
    {
        a[1] = -a[1];
        a[m2 + 1] = -a[m2 + 1];
        for (k = 1; k < m; k++)
        {
            for (j = 0; j < k; j++)
            {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            a[k1 + m2 + 1] = -a[k1 + m2 + 1];
        }
    }
}

void cftfsub(OsInt32 n,OsFloat *a,OsFloat *w)
{
    void cft1st(OsInt32 n, OsFloat *a, OsFloat *w);
    void cftmdl(OsInt32 n, OsInt32 l, OsFloat *a, OsFloat *w);
    OsInt32 j, j1, j2, j3, l;
    OsFloat x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    l = 2;
    if (n > 8)
    {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n)
        {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n)
    {
        for (j = 0; j < l; j += 2)
        {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
    }
    else
    {
        for (j = 0; j < l; j += 2)
        {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = a[j + 1] - a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] += a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}

void cftbsub(OsInt32 n, OsFloat *a, OsFloat *w)
{
    void cft1st(OsInt32 n, OsFloat *a, OsFloat *w);
    void cftmdl(OsInt32 n, OsInt32 l, OsFloat *a, OsFloat *w);
    OsInt32 j, j1, j2, j3, l;
    OsFloat x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    l = 2;
    if (n > 8)
    {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n)
        {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n)
    {
        for (j = 0; j < l; j += 2)
        {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = -a[j + 1] - a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = -a[j + 1] + a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i - x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i + x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i - x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i + x3r;
        }
    }
    else
    {
        for (j = 0; j < l; j += 2)
        {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = -a[j + 1] + a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] = -a[j + 1] - a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}

void cft1st(OsInt32 n, OsFloat *a, OsFloat *w)
{
    OsInt32 j, k1, k2;
    OsFloat wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    OsFloat x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    x0r = a[0] + a[2];
    x0i = a[1] + a[3];
    x1r = a[0] - a[2];
    x1i = a[1] - a[3];
    x2r = a[4] + a[6];
    x2i = a[5] + a[7];
    x3r = a[4] - a[6];
    x3i = a[5] - a[7];
    a[0] = x0r + x2r;
    a[1] = x0i + x2i;
    a[4] = x0r - x2r;
    a[5] = x0i - x2i;
    a[2] = x1r - x3i;
    a[3] = x1i + x3r;
    a[6] = x1r + x3i;
    a[7] = x1i - x3r;
    wk1r = w[2];
    x0r = a[8] + a[10];
    x0i = a[9] + a[11];
    x1r = a[8] - a[10];
    x1i = a[9] - a[11];
    x2r = a[12] + a[14];
    x2i = a[13] + a[15];
    x3r = a[12] - a[14];
    x3i = a[13] - a[15];
    a[8] = x0r + x2r;
    a[9] = x0i + x2i;
    a[12] = x2i - x0i;
    a[13] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[10] = wk1r * (x0r - x0i);
    a[11] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[14] = wk1r * (x0i - x0r);
    a[15] = wk1r * (x0i + x0r);
    k1 = 0;
    for (j = 16; j < n; j += 16)
    {
        k1 += 2;
        k2 = 2 * k1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        x0r = a[j] + a[j + 2];
        x0i = a[j + 1] + a[j + 3];
        x1r = a[j] - a[j + 2];
        x1i = a[j + 1] - a[j + 3];
        x2r = a[j + 4] + a[j + 6];
        x2i = a[j + 5] + a[j + 7];
        x3r = a[j + 4] - a[j + 6];
        x3i = a[j + 5] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 4] = wk2r * x0r - wk2i * x0i;
        a[j + 5] = wk2r * x0i + wk2i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 2] = wk1r * x0r - wk1i * x0i;
        a[j + 3] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 6] = wk3r * x0r - wk3i * x0i;
        a[j + 7] = wk3r * x0i + wk3i * x0r;
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        x0r = a[j + 8] + a[j + 10];
        x0i = a[j + 9] + a[j + 11];
        x1r = a[j + 8] - a[j + 10];
        x1i = a[j + 9] - a[j + 11];
        x2r = a[j + 12] + a[j + 14];
        x2i = a[j + 13] + a[j + 15];
        x3r = a[j + 12] - a[j + 14];
        x3i = a[j + 13] - a[j + 15];
        a[j + 8] = x0r + x2r;
        a[j + 9] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 12] = -wk2i * x0r - wk2r * x0i;
        a[j + 13] = -wk2i * x0i + wk2r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 10] = wk1r * x0r - wk1i * x0i;
        a[j + 11] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 14] = wk3r * x0r - wk3i * x0i;
        a[j + 15] = wk3r * x0i + wk3i * x0r;
    }
}

void cftmdl(OsInt32 n, OsInt32 l, OsFloat *a, OsFloat *w)
{
    OsInt32 j, j1, j2, j3, k, k1, k2, m, m2;
    OsFloat wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    OsFloat x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
    
    m = l << 2;
    for (j = 0; j < l; j += 2)
    {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x0r - x2r;
        a[j2 + 1] = x0i - x2i;
        a[j1] = x1r - x3i;
        a[j1 + 1] = x1i + x3r;
        a[j3] = x1r + x3i;
        a[j3 + 1] = x1i - x3r;
    }
    wk1r = w[2];
    for (j = m; j < l + m; j += 2)
    {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x2i - x0i;
        a[j2 + 1] = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j1] = wk1r * (x0r - x0i);
        a[j1 + 1] = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[j3] = wk1r * (x0i - x0r);
        a[j3 + 1] = wk1r * (x0i + x0r);
    }
    k1 = 0;
    m2 = 2 * m;
    for (k = m2; k < n; k += m2)
    {
        k1 += 2;
        k2 = 2 * k1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        for (j = k; j < l + k; j += 2)
        {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = wk2r * x0r - wk2i * x0i;
            a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        for (j = k + m; j < l + (k + m); j += 2)
        {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = -wk2i * x0r - wk2r * x0i;
            a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
    }
}

void rftfsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c)
{
    OsInt32 j, k, kk, ks, m;
    OsFloat wkr, wki, xr, xi, yr, yi;
    
    m = n >> 1;
    ks = 2 * nc / m;
    kk = 0;
    for (j = 2; j < m; j += 2)
    {
        k = n - j;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[j] - a[k];
        xi = a[j + 1] + a[k + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        a[j] -= yr;
        a[j + 1] -= yi;
        a[k] += yr;
        a[k + 1] -= yi;
    }
}

void rftbsub(OsInt32 n,OsFloat *a,OsInt32 nc,OsFloat *c)
{
    OsInt32 j, k, kk, ks, m;
    OsFloat wkr, wki, xr, xi, yr, yi;
    
    a[1] = -a[1];
    m = n >> 1;
    ks = 2 * nc / m;
    kk = 0;
    for (j = 2; j < m; j += 2)
    {
        k = n - j;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[j] - a[k];
        xi = a[j + 1] + a[k + 1];
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[j] -= yr;
        a[j + 1] = yi - a[j + 1];
        a[k] += yr;
        a[k + 1] = yi - a[k + 1];
    }
    a[m + 1] = -a[m + 1];
}


void dctsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c)
{
    OsInt32 j, k, kk, ks, m;
    OsFloat wkr, wki, xr;
    
    m = n >> 1;
    ks = nc / n;
    kk = 0;
    for (j = 1; j < m; j++)
    {
        k = n - j;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[j] - wkr * a[k];
        a[j] = wkr * a[j] + wki * a[k];
        a[k] = xr;
    }
    a[m] *= c[0];
}


void dstsub(OsInt32 n, OsFloat *a, OsInt32 nc, OsFloat *c)
{
    OsInt32 j, k, kk, ks, m;
    OsFloat wkr, wki, xr;
    
    m = n >> 1;
    ks = nc / n;
    kk = 0;
    for (j = 1; j < m; j++)
    {
        k = n - j;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[k] - wkr * a[j];
        a[k] = wkr * a[k] + wki * a[j];
        a[j] = xr;
    }
    a[m] *= c[0];
}

#define L_SUBFR 64
const OsFloat sincos_t[161] =
{
    0.0f,
    0.0245412290096282960f, 0.0490676760673522950f, 0.0735645666718482970f, 0.0980171412229537960f, 0.1224106773734092700f,
    0.1467304676771163900f, 0.1709618866443634000f, 0.1950903236865997300f, 0.2191012352705001800f, 0.2429801821708679200f,
    0.2667127549648284900f, 0.2902846634387970000f, 0.3136817514896392800f, 0.3368898630142211900f, 0.3598950505256652800f,
    0.3826834261417388900f, 0.4052413105964660600f, 0.4275550842285156300f, 0.4496113359928131100f, 0.4713967442512512200f,
    0.4928981959819793700f, 0.5141027569770813000f, 0.5349976420402526900f, 0.5555702447891235400f, 0.5758081674575805700f,
    0.5956993103027343800f, 0.6152315735816955600f, 0.6343932747840881300f, 0.6531728506088256800f, 0.6715589761734008800f,
    0.6895405650138855000f, 0.7071067690849304200f, 0.7242470979690551800f, 0.7409511208534240700f, 0.7572088241577148400f,
    0.7730104327201843300f, 0.7883464097976684600f, 0.8032075166702270500f, 0.8175848126411438000f, 0.8314695954322814900f,
    0.8448535799980163600f, 0.8577286005020141600f, 0.8700869679450988800f, 0.8819212913513183600f, 0.8932242989540100100f,
    0.9039893150329589800f, 0.9142097830772399900f, 0.9238795042037963900f, 0.9329928159713745100f, 0.9415440559387207000f,
    0.9495281577110290500f, 0.9569403529167175300f, 0.9637760519981384300f, 0.9700312614440918000f, 0.9757021069526672400f,
    0.9807852506637573200f, 0.9852776527404785200f, 0.9891765117645263700f, 0.9924795627593994100f, 0.9951847195625305200f,
    0.9972904324531555200f, 0.9987954497337341300f, 0.9996988177299499500f, 1.0000000000000000000f, 0.9996988177299499500f,
    0.9987954497337341300f, 0.9972904324531555200f, 0.9951847195625305200f, 0.9924795627593994100f, 0.9891765117645263700f,
    0.9852776527404785200f, 0.9807852506637573200f, 0.9757021069526672400f, 0.9700312614440918000f, 0.9637760519981384300f,
    0.9569403529167175300f, 0.9495281577110290500f, 0.9415440559387207000f, 0.9329928159713745100f, 0.9238795042037963900f,
    0.9142097830772399900f, 0.9039893150329589800f, 0.8932242989540100100f, 0.8819212913513183600f, 0.8700869679450988800f,
    0.8577286005020141600f, 0.8448535799980163600f, 0.8314695954322814900f, 0.8175848126411438000f, 0.8032075166702270500f,
    0.7883464097976684600f, 0.7730104327201843300f, 0.7572088241577148400f, 0.7409511208534240700f, 0.7242470979690551800f,
    0.7071067690849304200f, 0.6895405650138855000f, 0.6715589761734008800f, 0.6531728506088256800f, 0.6343932747840881300f,
    0.6152315735816955600f, 0.5956993103027343800f, 0.5758081674575805700f, 0.5555702447891235400f, 0.5349976420402526900f,
    0.5141027569770813000f, 0.4928981959819793700f, 0.4713967442512512200f, 0.4496113359928131100f, 0.4275550842285156300f,
    0.4052413105964660600f, 0.3826834261417388900f, 0.3598950505256652800f, 0.3368898630142211900f, 0.3136817514896392800f,
    0.2902846634387970000f, 0.2667127549648284900f, 0.2429801821708679200f, 0.2191012352705001800f, 0.1950903236865997300f,
    0.1709618866443634000f, 0.1467304676771163900f, 0.1224106773734092700f, 0.0980171412229537960f, 0.0735645666718482970f,
    0.0490676723420619960f, 0.0245412290096282960f, -0.00000000041020686847303978f, -0.0245412290096282960f, -0.0490676760673522950f,
    -0.0735645666718482970f, -0.0980171412229537960f, -0.1224106773734092700f, -0.1467304676771163900f, -0.1709618866443634000f,
    -0.1950903236865997300f, -0.2191012352705001800f, -0.2429801821708679200f, -0.2667127549648284900f, -0.2902846634387970000f,
    -0.3136817514896392800f, -0.3368898630142211900f, -0.3598950505256652800f, -0.3826834261417388900f, -0.4052413105964660600f,
    -0.4275550842285156300f, -0.4496113359928131100f, -0.4713967442512512200f, -0.4928981959819793700f, -0.5141027569770813000f,
    -0.5349976420402526900f, -0.5555702447891235400f, -0.5758081674575805700f, -0.5956993103027343800f, -0.6152315735816955600f,
    -0.6343932747840881300f, -0.6531728506088256800f, -0.6715589761734008800f, -0.6895405650138855000f, -0.7071067690849304200f
};

#define N_MAX       256

/*---------------------------------------------------------------------*
 *  fft_rel()
 *
 *  Computes the split-radix FFT in place for the real-valued
 *  signal x of length n.  The algorithm has been ported from
 *  the Fortran code of [1].
 *
 *  The function  needs sine and cosine tables t_sin and t_cos,
 *  and the constant N_MAX.  The table  entries  are defined as
 *  sin(2*pi*i) and cos(2*pi*i) for i = 0, 1, ..., N_MAX-1. The
 *  implementation  assumes  that any entry  will not be needed
 *  outside the tables. Therefore, N_MAX and n must be properly
 *  set.  The function has been tested  with the values n = 16,
 *  32, 64, 128, 256, and N_MAX = 1280.
 *
 *  References
 *  [1] H.V. Sorensen,  D.L. Jones, M.T. Heideman, C.S. Burrus,
 *      "Real-valued fast  Fourier transform  algorithm,"  IEEE
 *      Trans. on Signal Processing,  Vol.35, No.6, pp 849-863,
 *      1987.
 *
 *  OUTPUT
 *      x[0:n-1]  Transform coeffients in the order re[0], re[1],
 *                ..., re[n/2], im[n/2-1], ..., im[1].
 *---------------------------------------------------------------------*/

/* i/o: input/output vector    */
/* i  : vector length          */
/* i  : log2 of vector length  */
void fft_rel(OsFloat x[],const short n,const short m)
{
    short  i, j, k, n1, n2, n4;
    short  step;
    OsFloat  xt, t1, t2;
    OsFloat *x0, *x1, *x2, *s, *c;
    OsFloat *xi2, *xi3, *xi4, *xi1;

    /*-----------------------------------------------------------------*
    * Digit reverse counter
    *-----------------------------------------------------------------*/

    j = 0;
    x0 = &x[0];
    for (i = 0; i < n-1; i++)
    {
        if (i < j)
        {
            xt   = x[j];
            x[j] = *x0;
            *x0  = xt;
        }
        x0++;
        k = n/2;
        while (k <= j)
        {
            j -= k;
            k  = k>>1;
        }
        j += k;
    }

    /*-----------------------------------------------------------------*
    * Length two butterflies
    *-----------------------------------------------------------------*/

    x0 = &x[0];
    x1 = &x[1];
    for (i = 0; i < n/2; i++)
    {
        xt  = *x0;
        *x0 = xt + *x1;
        *x1 = xt - *x1;
        x0++; x0++;
        x1++; x1++;
    }

    /*-----------------------------------------------------------------*
    * Other butterflies
    *
    * The implementation described in [1] has been changed by using
    * table lookup for evaluating sine and cosine functions.  The
    * variable ind and its increment step are needed to access table
    * entries.  Note that this implementation assumes n4 to be so
    * small that ind will never exceed the table.  Thus the input
    * argument n and the constant N_MAX must be set properly.
    *-----------------------------------------------------------------*/

    n2 = 1;
    for (k = 2; k <= m; k++)
    {
        n4 = n2;
        n2 = n4<<1;
        n1 = n2<<1;
        step = N_MAX/n1;

        x0 = x;
        x1 = x + n2;
        x2 = x + n2 + n4;
        for (i = 0; i < n; i += n1)
        {
            xt = *x0;               /* xt = x[i];   */
            *x0 =  xt + *x1;        /* x[i] = xt + x[i+n2];    */
            *x1  = xt - *x1;        /* x[i+n2] = xt - x[i+n2];      */
            *x2 = -*x2;             /* x[i+n2+n4] = -x[i+n2+n4];     */

            s = (OsFloat *)sincos_t + step;
            c = s + 64;
            xi1 = x + i + 1;
            xi3 = xi1 + n2 ;
            xi2 = xi3 - 2;
            xi4 = xi1 + n1 - 2;

            for (j = 1; j < n4; j++)
            {
                t1  =  *xi3**c + *xi4**s;    /* t1 = *xi3**(pt_c+ind) + *xi4**(pt_s+ind);   */
                t2 = *xi3**s - *xi4**c;      /* t2 = *xi3**(pt_s+ind) - *xi4**(pt_c+ind);     */
                *xi4 = *xi2 - t2;
                *xi3 = -*xi2 - t2;
                *xi2 =  *xi1 - t1;
                *xi1 =  *xi1 + t1;

                xi4--; xi2--;
                xi3++; xi1++;
                c += step; s+= step;          /* autoincrement by ar0 */
            }
            x0 += n1;
            x1 += n1;
            x2 += n1;
        }
    }
    return;
}

#define INV_SQR2  0.70710676908493f

/*---------------------------------------------------------------------*
* ifft_rel()
*
* Calculate the inverse FFT of a real signal
*
* Based on the FORTRAN code from the article "Real-valued Fast Fourier Transform Algorithms"
* by Sorensen, ... in IEEE Trans. on ASSP, Vol. ASSP-35, No. June 6th 1987.
*
* Input: the io[] signal containing the spectrum in the following order :
*
* Re[0], Re[1], ..  Re[n/2], Im[n/2-1], .. Im[1]
*---------------------------------------------------------------------*/

/* i/o: input/output vector   */
/* i  : vector length         */
/* i  : log2 of vector length */
void ifft_rel(OsFloat io[],const short n,const short m)
{
    short i, j, k;
    short step;
    short n2, n4, n8, i0;
    short is, id;
    OsFloat *x,*xi0, *xi1, *xi2, *xi3, *xi4, *xup1, *xdn6, *xup3, *xdn8;
    OsFloat xt;
    OsFloat r1;
    OsFloat t1, t2, t3, t4, t5;
    OsFloat *s, *c, *s3, *c3;

    OsFloat cc1, cc3, ss1, ss3;


    /*-----------------------------------------------------------------*
    * ifft
    *-----------------------------------------------------------------*/

    x     = &io[-1];
    n2 = 2*n;
    for (k=1; k<m; k++)
    {
        is = 0;
        id = n2;
        n2 = n2 >> 1;
        n4 = n2 >> 2;
        n8 = n4 >> 1;
        while (is < n-1)
        {
            xi1 = x + is + 1;
            xi2 = xi1 + n4;
            xi3 = xi2 + n4;
            xi4 = xi3 + n4;

            for (i=is; i<n; i+= id)
            {
                t1 = *xi1 - *xi3;
                *xi1 += *xi3;
                *xi2 = 2.0f**xi2;
                *xi3 = t1- 2.0f**xi4;
                *xi4 = t1 + 2.0f**xi4;
                if (n4 != 1)
                {
                    t1 = (*(xi2+n8) - *(xi1+n8))*INV_SQR2;
                    t2 = (*(xi4+n8) + *(xi3+n8))*INV_SQR2;

                    *(xi1+n8) += *(xi2+n8);
                    *(xi2+n8) = *(xi4+n8) - *(xi3+n8);
                    *(xi3+n8) = (OsFloat)(2.0f * (-t2-t1));
                    *(xi4+n8) = (OsFloat)(2.0f * (-t2+t1));
                }
                xi1 += id;
                xi2 += id;
                xi3 += id;
                xi4 += id;
            }
            is = 2*id - n2;
            id  = 4*id;
        }
        step  = N_MAX/n2;

        s = (OsFloat *)sincos_t + step;
        c = s + 64;
        s3 = (OsFloat *)sincos_t + 3*step;
        c3 = s3 + 64;
        for (j=2; j<=n8; j++)
        {
            cc1 = *c ;
            ss1 = *s;
            cc3 = *c3;
            ss3 = *s3;

            is  = 0;
            id  = 2 * n2;

            c += step;
            s += step;

            c3 += 3*step;
            s3 += 3*step;
            while (is < n-1)
            {
                xup1 = x + j + is;
                xup3 = xup1 + 2*n4;
                xdn6 = xup3 - 2*j +2;
                xdn8 = xdn6 + 2*n4;

                for (i=is; i<n; i+=id)
                {
                    t1 = *xup1 - *xdn6;
                    *xup1 = *xup1 + *xdn6;
                    xup1 += n4;
                    xdn6 -= n4;

                    t2 = *xdn6 - *xup1;
                    *xdn6 = *xup1 + *xdn6;

                    xdn6 += n4;
                    t3 = *xdn8 + *xup3;
                    *xdn6  = *xdn8 - *xup3;

                    xup3 += n4;
                    xdn8 -= n4;

                    t4 =  *xup3 + *xdn8;
                    *xup1= *xup3 - *xdn8;

                    t5 = t1 - t4;
                    t1 = t1 + t4;
                    t4 = t2 - t3;
                    t2 = t2 + t3;
                    *xup3 = t1*cc3 - t2*ss3;
                    xup3 -= n4;
                    *xup3 = t5*cc1 + t4*ss1;
                    *xdn8 = -t4*cc1 + t5*ss1;

                    xdn8 += n4;
                    *xdn8 = t2*cc3 + t1*ss3;

                    xup1 -= n4;
                    xup1 += id;
                    xup3 += id;
                    xdn6 += id;
                    xdn8 += id;
                }
                is = 2*id - n2;
                id = 4*id;
            }
        }
    }

    /*-----------------------------------------------------------------*
    * Length two butterflies
    *-----------------------------------------------------------------*/

    is = 1;
    id = 4;
    while (is < n)
    {
        xi0 = x + is ;
        xi1 = xi0 + 1;

        for (i0=is; i0<=n; i0+=id)
        {
            r1 = *xi0;
            *xi0= r1 + *xi1;
            *xi1 = r1 - *xi1;
            xi0 += id;
            xi1 += id;
        }
        is = 2*id - 1;
        id = 4*id;
    }

    /*-----------------------------------------------------------------*
    * Digit reverse counter
    *-----------------------------------------------------------------*/

    j = 1;
    for (i=1; i<n; i++)
    {
        if (i < j)
        {
            xt = x[j];
            x[j] = x[i];
            x[i] = xt;
        }
        k = n >> 1;
        while (k < j)
        {
            j = j - k;
            k = k >> 1;
        }
        j = j + k;
    }
    /*-----------------------------------------------------------------*
    * Normalization
    *-----------------------------------------------------------------*/
    for (i=1; i<=n; i++)
    {
        x[i] = x[i] / (OsFloat) n;
    }
    return;
}


/*-------------------------------------------------------------------*
* DCTII, length of 64 is supported only
*-------------------------------------------------------------------*/
/* i  : time domain input       */
/* o  : transform domain output */
void dct2(const OsFloat in[],OsFloat out[])
{
    short i;
    OsFloat x[L_SUBFR];
    const OsFloat *weight_dct;

    /* Re-order elements */
    for(i = 0; i < L_SUBFR/2; i++)
    {
        x[i] = in[i*2];
        x[L_SUBFR-1-i] = in[i*2+1];
    }

    /* set pointer of DFT weights to the cos() values */
    weight_dct = sincos_t + 64;

    /* compute DCT via FFT */
    fft_rel(x,L_SUBFR,6);

    /* multiply weights (complex values) with fft coefficients (complex values) */
    /* re-order and prepare the output */
    out[0] = x[0];
    out[L_SUBFR/2] = weight_dct[L_SUBFR/2] * x[L_SUBFR/2];  /* Im(x[L_SUBFR/2]) is always 0 */

    for(i = 1; i < L_SUBFR/2; i++)
    {
        out[i] = weight_dct[i]*x[i] + weight_dct[L_SUBFR-i]*x[L_SUBFR-i];
        out[i+L_SUBFR/2] = weight_dct[i+L_SUBFR/2]*x[L_SUBFR/2-i] - weight_dct[L_SUBFR/2-i]*x[L_SUBFR/2+i];
    }
    return;
}

/*-------------------------------------------------------------------*
* iDCTII, length of 64 is supported only
* in  : transform domain input
* out : time domain output
*-------------------------------------------------------------------*/
void idct2(const OsFloat in[],OsFloat out[])
{
    short i;
    OsFloat x[2*L_SUBFR];
    const OsFloat *weight_dct;

    /* set pointer of DFT weights to the cos() values */
    weight_dct = sincos_t + 64;

    /* multiply weights (complex values) with DCT coefficients (real values) */
    x[0] = in[0];
    x[L_SUBFR] = -1*x[0];

    for(i = 1; i < L_SUBFR; i++)
    {
        x[i] = 2 * weight_dct[i] * in[i];
        x[2*L_SUBFR-i] = 2 * weight_dct[L_SUBFR-i] * in[i];
    }

    /* compute iDCT via iFFT */
    ifft_rel(x,2*L_SUBFR, 6+1);

    /* re-order elements */
    for(i = 0; i < L_SUBFR/2; i++)
    {
        out[i*2+1] = x[i*2+1];
        out[i*2]   = x[2*L_SUBFR-1-i*2];
    }
    return;
}

/*************************************** speex_dft ****************************************/

struct kft_t
{
    OsInt32 n;
    OsFloat *trigcache;
    OsInt32 *splitcache;
    OsFloat *x;
};

extern void kiss_drft_forward(kft_t *l,OsFloat *data);
extern void kiss_drft_backward(kft_t *l,OsFloat *data);
extern void kiss_drft_init(kft_t *l,OsInt32 n);
extern void kiss_drft_clear(kft_t *l);

kft_t* kiss_fft_init(OsInt32 size)
{
    kft_t *fft = (kft_t*)icm_alloc(sizeof(kft_t),1);
    kiss_drft_init(fft,size);
    fft->x = (OsFloat*)icm_alloc(size,sizeof(OsFloat));
    return fft;
}

void kiss_fft_uninit(kft_t *fft)
{
    kiss_drft_clear(fft);
    icm_free(fft->x);
    icm_free(fft);
}

void kiss_fft(kft_t *fft,OsFloat *x,complex *fx)
{
    OsInt32 bins = fft->n/2+1;
    OsInt32 i,l = 0;

    memcpy(fft->x,x,fft->n*sizeof(OsFloat));
    kiss_drft_forward(fft,fft->x);
    fx[0].re = fft->x[0]; fx[0].im = 0;
    if(0 == (fft->n % 2))
    {
        for(i = 1,l = 1; i < bins-1; i++,l += 2)
        {
            fx[i].re = fft->x[l];
            fx[i].im = fft->x[l+1];
        }
        fx[i].re = fft->x[l]; fx[i].im = 0;
        //for(i = bins, l = bins-2; i < n; i++,l--)
        //{
        //    fx[i].re = fx[l].re;
        //    fx[i].im = -fx[l].im;
        //}
    }
    else
    {
        for(i = 1,l = 1; i < bins; i++,l += 2)
        {
            fx[i].re = fft->x[l];
            fx[i].im = fft->x[l+1];
        }
        //for(i = bins, l = bins-1; i < n; i++,l--)
        //{
        //    fx[i].re = fx[l].re;
        //    fx[i].im = -fx[l].im;
        //}
    }
}

void kiss_ifft(kft_t *fft,complex *fx,OsFloat *x)
{
    OsInt32 bins = fft->n/2+1;
    OsFloat scale = 1.0f / fft->n;
    OsInt32 i,l;

    // 转为rdft需要的排列，同时做一次scale
    x[0] = fx[0].re*scale;
    if(0 == fft->n%2)
    {
        for(i = 1,l = 1; i < bins-1; i++,l += 2)
        {
            x[l] = fx[i].re*scale;
            x[l+1] = fx[i].im*scale;
        }
        x[l] = fx[i].re*scale;
    }
    else
    {
        for(i = 1,l = 1; i < bins; i++,l += 2)
        {
            x[l] = fx[i].re*scale;
            x[l+1] = fx[i].im*scale;
        }
    }
    kiss_drft_backward(fft,x);
}

void kiss_fft2(kft_t *fft,OsFloat *x)
{
    kiss_drft_forward(fft,x);
}

void kiss_ifft2(kft_t *fft,OsFloat *x)
{
    OsInt32 n = fft->n;
    OsFloat scale = 1.0f / n;
    OsInt32 i;

    for (i = 0; i < n; i++)
    {
        x[i] *= scale;
    }
    kiss_drft_backward(fft,x);
}

static void drfti1(OsInt32 n,OsFloat *wa,OsInt32 *ifac)
{
    static OsInt32 ntryh[4] = { 4,2,3,5 };
    static OsFloat tpi = 6.28318530717958648f;
    OsFloat arg,argh,argld,fi;
    OsInt32 ntry=0,i,j=-1;
    OsInt32 k1, l1, l2, ib;
    OsInt32 ld, ii, ip, is, nq, nr;
    OsInt32 ido, ipm, nfm1;
    OsInt32 nl=n;
    OsInt32 nf=0;
 L101:
    j++;
    if (j < 4)
        ntry=ntryh[j];
    else
        ntry+=2;
 L104:
    nq=nl/ntry;
    nr=nl-ntry*nq;
    if (nr!=0) goto L101;
    nf++;
    ifac[nf+1]=ntry;
    nl=nq;
    if(ntry!=2)goto L107;
    if(nf==1)goto L107;
    for (i=1;i<nf;i++)
    {
        ib=nf-i+1;
        ifac[ib+1]=ifac[ib];
    }
    ifac[2] = 2;
 L107:
    if(nl!=1)goto L104;
    ifac[0]=n;
    ifac[1]=nf;
    argh=tpi/n;
    is=0;
    nfm1=nf-1;
    l1=1;
    if(nfm1==0)return;
    for (k1=0;k1<nfm1;k1++)
    {
        ip=ifac[k1+2];
        ld=0;
        l2=l1*ip;
        ido=n/l2;
        ipm=ip-1;
        for (j=0;j<ipm;j++)
        {
            ld += l1;
            i = is;
            argld = (OsFloat)ld*argh;
            fi = 0.f;
            for (ii=2;ii<ido;ii+=2)
            {
	            fi += 1.f;
	            arg = fi*argld;
	            wa[i++] = cosf(arg);
	            wa[i++] = sinf(arg);
            }
            is+=ido;
        }
        l1=l2;
    }
}

static void fdrffti(OsInt32 n,OsFloat *wsave,OsInt32 *ifac)
{
    if (n == 1) return;
    drfti1(n, wsave+n, ifac);
}

static void dradf2(OsInt32 ido,OsInt32 l1,OsFloat *cc,OsFloat *ch,OsFloat *wa1)
{
    OsInt32 i,k;
    OsFloat ti2,tr2;
    OsInt32 t0,t1,t2,t3,t4,t5,t6;

    t1=0;
    t0=(t2=l1*ido);
    t3=ido<<1;
    for(k=0;k<l1;k++)
    {
        ch[t1<<1]=cc[t1]+cc[t2];
        ch[(t1<<1)+t3-1]=cc[t1]-cc[t2];
        t1+=ido;
        t2+=ido;
    }

    if(ido<2)return;
    if(ido==2)goto L105;

    t1=0;
    t2=t0;
    for(k=0;k<l1;k++)
    {
        t3=t2;
        t4=(t1<<1)+(ido<<1);
        t5=t1;
        t6=t1+t1;
        for(i=2;i<ido;i+=2)
        {
            t3+=2;
            t4-=2;
            t5+=2;
            t6+=2;
            tr2=wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
            ti2=wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
            ch[t6]=cc[t5]+ti2;
            ch[t4]=ti2-cc[t5];
            ch[t6-1]=cc[t5-1]+tr2;
            ch[t4-1]=cc[t5-1]-tr2;
        }
        t1+=ido;
        t2+=ido;
    }
    if(ido%2==1)return;
 L105:
    t3=(t2=(t1=ido)-1);
    t2+=t0;
    for(k=0;k<l1;k++)
    {
        ch[t1]=-cc[t2];
        ch[t1-1]=cc[t3];
        t1+=ido<<1;
        t2+=ido;
        t3+=ido;
    }
}

static void dradf4(OsInt32 ido,OsInt32 l1,OsFloat *cc,OsFloat *ch,OsFloat *wa1,OsFloat *wa2,OsFloat *wa3)
{
    static OsFloat hsqt2 = .70710678118654752f;
    OsInt32 i,k,t0,t1,t2,t3,t4,t5,t6;
    OsFloat ci2,ci3,ci4,cr2,cr3,cr4,ti1,ti2,ti3,ti4,tr1,tr2,tr3,tr4;
    t0=l1*ido;
  
    t1=t0;
    t4=t1<<1;
    t2=t1+(t1<<1);
    t3=0;

    for(k=0;k<l1;k++)
    {
        tr1=cc[t1]+cc[t2];
        tr2=cc[t3]+cc[t4];

        ch[t5=t3<<2]=tr1+tr2;
        ch[(ido<<2)+t5-1]=tr2-tr1;
        ch[(t5+=(ido<<1))-1]=cc[t3]-cc[t4];
        ch[t5]=cc[t2]-cc[t1];

        t1+=ido;
        t2+=ido;
        t3+=ido;
        t4+=ido;
    }

    if(ido<2)return;
    if(ido==2)goto L105;

    t1=0;
    for(k=0;k<l1;k++)
    {
        t2=t1;
        t4=t1<<2;
        t5=(t6=ido<<1)+t4;
        for(i=2;i<ido;i+=2)
        {
            t3=(t2+=2);
            t4+=2;
            t5-=2;

            t3+=t0;
            cr2=wa1[i-2]*cc[t3-1]+wa1[i-1]*cc[t3];
            ci2=wa1[i-2]*cc[t3]-wa1[i-1]*cc[t3-1];
            t3+=t0;
            cr3=wa2[i-2]*cc[t3-1]+wa2[i-1]*cc[t3];
            ci3=wa2[i-2]*cc[t3]-wa2[i-1]*cc[t3-1];
            t3+=t0;
            cr4=wa3[i-2]*cc[t3-1]+wa3[i-1]*cc[t3];
            ci4=wa3[i-2]*cc[t3]-wa3[i-1]*cc[t3-1];

            tr1=cr2+cr4;
            tr4=cr4-cr2;
            ti1=ci2+ci4;
            ti4=ci2-ci4;

            ti2=cc[t2]+ci3;
            ti3=cc[t2]-ci3;
            tr2=cc[t2-1]+cr3;
            tr3=cc[t2-1]-cr3;

            ch[t4-1]=tr1+tr2;
            ch[t4]=ti1+ti2;

            ch[t5-1]=tr3-ti4;
            ch[t5]=tr4-ti3;

            ch[t4+t6-1]=ti4+tr3;
            ch[t4+t6]=tr4+ti3;

            ch[t5+t6-1]=tr2-tr1;
            ch[t5+t6]=ti1-ti2;
        }
        t1+=ido;
    }
    if(ido&1)return;
 L105:
    t2=(t1=t0+ido-1)+(t0<<1);
    t3=ido<<2;
    t4=ido;
    t5=ido<<1;
    t6=ido;

    for(k=0;k<l1;k++)
    {
        ti1=-hsqt2*(cc[t1]+cc[t2]);
        tr1=hsqt2*(cc[t1]-cc[t2]);

        ch[t4-1]=tr1+cc[t6-1];
        ch[t4+t5-1]=cc[t6-1]-tr1;

        ch[t4]=ti1-cc[t1+t0];
        ch[t4+t5]=ti1+cc[t1+t0];

        t1+=ido;
        t2+=ido;
        t4+=t3;
        t6+=ido;
    }
}

static void dradfg(OsInt32 ido,OsInt32 ip,OsInt32 l1,OsInt32 idl1,OsFloat *cc,OsFloat *c1,OsFloat *c2,OsFloat *ch,OsFloat *ch2,OsFloat *wa)
{
    static OsFloat tpi=6.283185307179586f;
    OsInt32 idij,ipph,i,j,k,l,ic,ik,is;
    OsInt32 t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    OsFloat dc2,ai1,ai2,ar1,ar2,ds2;
    OsInt32 nbd;
    OsFloat dcp,arg,dsp,ar1h,ar2h;
    OsInt32 idp2,ipp2;

    arg =tpi/(OsFloat)ip;
    dcp = cosf(arg);
    dsp = sinf(arg);
    ipph = (ip+1)>>1;
    ipp2 = ip;
    idp2 = ido;
    nbd = (ido-1)>>1;
    t0  = l1*ido;
    t10 = ip*ido;

    if(ido == 1) goto L119;
    for(ik = 0; ik < idl1; ik++) ch2[ik]=c2[ik];

    t1=0;
    for(j=1;j<ip;j++)
    {
        t1+=t0;
        t2=t1;
        for(k=0;k<l1;k++)
        {
            ch[t2]=c1[t2];
            t2+=ido;
        }
    }
    is=-ido;
    t1=0;
    if(nbd>l1)
    {
        for(j=1;j<ip;j++)
        {
            t1+=t0;
            is+=ido;
            t2= -ido+t1;
            for(k=0;k<l1;k++)
            {
                idij=is-1;
                t2+=ido;
                t3=t2;
                for(i=2;i<ido;i+=2)
                {
                    idij+=2;
                    t3+=2;
                    ch[t3-1]=wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
                    ch[t3]=wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
                }
            }
        }
    }
    else
    {
        for(j=1;j<ip;j++)
        {
            is+=ido;
            idij=is-1;
            t1+=t0;
            t2=t1;
            for(i=2;i<ido;i+=2)
            {
                idij+=2;
                t2+=2;
                t3=t2;
                for(k=0;k<l1;k++)
                {
                    ch[t3-1]=wa[idij-1]*c1[t3-1]+wa[idij]*c1[t3];
                    ch[t3]=wa[idij-1]*c1[t3]-wa[idij]*c1[t3-1];
                    t3+=ido;
                }
            }
        }
    }
    t1=0;
    t2=ipp2*t0;
    if(nbd<l1)
    {
        for(j=1;j<ipph;j++)
        {
            t1+=t0;
            t2-=t0;
            t3=t1;
            t4=t2;
            for(i=2;i<ido;i+=2)
            {
                t3+=2;
                t4+=2;
                t5=t3-ido;
                t6=t4-ido;
                for(k=0;k<l1;k++)
                {
                    t5+=ido;
                    t6+=ido;
                    c1[t5-1]=ch[t5-1]+ch[t6-1];
                    c1[t6-1]=ch[t5]-ch[t6];
                    c1[t5]=ch[t5]+ch[t6];
                    c1[t6]=ch[t6-1]-ch[t5-1];
                }
            }
        }
    }
    else
    {
        for(j=1;j<ipph;j++)
        {
            t1+=t0;
            t2-=t0;
            t3=t1;
            t4=t2;
            for(k=0;k<l1;k++)
            {
                t5=t3;
                t6=t4;
                for(i=2;i<ido;i+=2)
                {
                    t5+=2;
                    t6+=2;
                    c1[t5-1]=ch[t5-1]+ch[t6-1];
                    c1[t6-1]=ch[t5]-ch[t6];
                    c1[t5]=ch[t5]+ch[t6];
                    c1[t6]=ch[t6-1]-ch[t5-1];
                }
                t3+=ido;
                t4+=ido;
            }
        }
    }

L119:
    for(ik=0;ik<idl1;ik++) c2[ik]=ch2[ik];

    t1=0;
    t2=ipp2*idl1;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1-ido;
        t4=t2-ido;
        for(k=0;k<l1;k++)
        {
            t3+=ido;
            t4+=ido;
            c1[t3]=ch[t3]+ch[t4];
            c1[t4]=ch[t4]-ch[t3];
        }
    }

    ar1=1.f;
    ai1=0.f;
    t1=0;
    t2=ipp2*idl1;
    t3=(ip-1)*idl1;
    for(l=1;l<ipph;l++)
    {
        t1+=idl1;
        t2-=idl1;
        ar1h=dcp*ar1-dsp*ai1;
        ai1=dcp*ai1+dsp*ar1;
        ar1=ar1h;
        t4=t1;
        t5=t2;
        t6=t3;
        t7=idl1;

        for(ik=0;ik<idl1;ik++)
        {
            ch2[t4++]=c2[ik]+ar1*c2[t7++];
            ch2[t5++]=ai1*c2[t6++];
        }

        dc2=ar1;
        ds2=ai1;
        ar2=ar1;
        ai2=ai1;

        t4=idl1;
        t5=(ipp2-1)*idl1;
        for(j=2;j<ipph;j++)
        {
            t4+=idl1;
            t5-=idl1;

            ar2h=dc2*ar2-ds2*ai2;
            ai2=dc2*ai2+ds2*ar2;
            ar2=ar2h;

            t6=t1;
            t7=t2;
            t8=t4;
            t9=t5;
            for(ik=0;ik<idl1;ik++)
            {
                ch2[t6++]+=ar2*c2[t8++];
                ch2[t7++]+=ai2*c2[t9++];
            }
        }
    }

    t1=0;
    for(j=1;j<ipph;j++)
    {
        t1+=idl1;
        t2=t1;
        for(ik=0;ik<idl1;ik++)ch2[ik]+=c2[t2++];
    }

    if(ido<l1)goto L132;

    t1=0;
    t2=0;
    for(k=0;k<l1;k++)
    {
        t3=t1;
        t4=t2;
        for(i=0;i<ido;i++)cc[t4++]=ch[t3++];
        t1+=ido;
        t2+=t10;
    }
    goto L135;
L132:
    for(i=0;i<ido;i++)
    {
        t1=i;
        t2=i;
        for(k=0;k<l1;k++)
        {
            cc[t2]=ch[t1];
            t1+=ido;
            t2+=t10;
        }
    }

L135:
    t1=0;
    t2=ido<<1;
    t3=0;
    t4=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t2;
        t3+=t0;
        t4-=t0;

        t5=t1;
        t6=t3;
        t7=t4;

        for(k=0;k<l1;k++)
        {
            cc[t5-1]=ch[t6];
            cc[t5]=ch[t7];
            t5+=t10;
            t6+=ido;
            t7+=ido;
        }
    }

    if(ido==1)return;
    if(nbd<l1)goto L141;

    t1=-ido;
    t3=0;
    t4=0;
    t5=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t2;
        t3+=t2;
        t4+=t0;
        t5-=t0;
        t6=t1;
        t7=t3;
        t8=t4;
        t9=t5;
        for(k=0;k<l1;k++)
        {
            for(i=2;i<ido;i+=2)
            {
                ic=idp2-i;
                cc[i+t7-1]=ch[i+t8-1]+ch[i+t9-1];
                cc[ic+t6-1]=ch[i+t8-1]-ch[i+t9-1];
                cc[i+t7]=ch[i+t8]+ch[i+t9];
                cc[ic+t6]=ch[i+t9]-ch[i+t8];
            }
            t6+=t10;
            t7+=t10;
            t8+=ido;
            t9+=ido;
        }
    }
    return;

L141:
    t1=-ido;
    t3=0;
    t4=0;
    t5=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t2;
        t3+=t2;
        t4+=t0;
        t5-=t0;
        for(i=2;i<ido;i+=2)
        {
            t6=idp2+t1-i;
            t7=i+t3;
            t8=i+t4;
            t9=i+t5;
            for(k=0;k<l1;k++)
            {
                cc[t7-1]=ch[t8-1]+ch[t9-1];
                cc[t6-1]=ch[t8-1]-ch[t9-1];
                cc[t7]=ch[t8]+ch[t9];
                cc[t6]=ch[t9]-ch[t8];
                t6+=t10;
                t7+=t10;
                t8+=ido;
                t9+=ido;
            }
        }
    }
}

static void drftf1(OsInt32 n,OsFloat *c,OsFloat *ch,OsFloat *wa,OsInt32 *ifac)
{
    OsInt32 i,k1,l1,l2;
    OsInt32 na,kh,nf;
    OsInt32 ip,iw,ido,idl1,ix2,ix3;

    nf=ifac[1];
    na=1;
    l2=n;
    iw=n;

    for(k1=0;k1<nf;k1++)
    {
        kh=nf-k1;
        ip=ifac[kh+1];
        l1=l2/ip;
        ido=n/l2;
        idl1=ido*l1;
        iw-=(ip-1)*ido;
        na=1-na;

        if(ip!=4)goto L102;

        ix2=iw+ido;
        ix3=ix2+ido;
        if(na!=0)
            dradf4(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1);
        else
            dradf4(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1);
        goto L110;
L102:
        if(ip!=2)goto L104;
        if(na!=0)goto L103;

        dradf2(ido,l1,c,ch,wa+iw-1);
        goto L110;
L103:
        dradf2(ido,l1,ch,c,wa+iw-1);
        goto L110;
L104:
        if(ido==1)na=1-na;
        if(na!=0)goto L109;

        dradfg(ido,ip,l1,idl1,c,c,c,ch,ch,wa+iw-1);
        na=1;
        goto L110;
L109:
        dradfg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa+iw-1);
        na=0;
L110:
        l2=l1;
    }
    if(na == 1) return;
    for(i = 0; i < n; i++) c[i] = ch[i];
}

static void dradb2(OsInt32 ido,OsInt32 l1,OsFloat *cc,OsFloat *ch,OsFloat *wa1)
{
    OsInt32 i,k,t0,t1,t2,t3,t4,t5,t6;
    OsFloat ti2,tr2;

    t0=l1*ido;

    t1=0;
    t2=0;
    t3=(ido<<1)-1;
    for(k=0;k<l1;k++)
    {
        ch[t1]=cc[t2]+cc[t3+t2];
        ch[t1+t0]=cc[t2]-cc[t3+t2];
        t2=(t1+=ido)<<1;
    }

    if(ido<2)return;
    if(ido==2)goto L105;

    t1=0;
    t2=0;
    for(k=0;k<l1;k++)
    {
        t3=t1;
        t5=(t4=t2)+(ido<<1);
        t6=t0+t1;
        for(i=2;i<ido;i+=2)
        {
            t3+=2;
            t4+=2;
            t5-=2;
            t6+=2;
            ch[t3-1]=cc[t4-1]+cc[t5-1];
            tr2=cc[t4-1]-cc[t5-1];
            ch[t3]=cc[t4]-cc[t5];
            ti2=cc[t4]+cc[t5];
            ch[t6-1]=wa1[i-2]*tr2-wa1[i-1]*ti2;
            ch[t6]=wa1[i-2]*ti2+wa1[i-1]*tr2;
        }
        t2=(t1+=ido)<<1;
    }

    if(ido%2==1)return;

L105:
    t1=ido-1;
    t2=ido-1;
    for(k=0;k<l1;k++)
    {
        ch[t1]=cc[t2]+cc[t2];
        ch[t1+t0]=-(cc[t2+1]+cc[t2+1]);
        t1+=ido;
        t2+=ido<<1;
    }
}

static void dradb3(OsInt32 ido,OsInt32 l1,OsFloat *cc,OsFloat *ch,OsFloat *wa1,OsFloat *wa2)
{
    static OsFloat taur = -.5f;
    static OsFloat taui = .8660254037844386f;
    OsInt32 i,k,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    OsFloat ci2,ci3,di2,di3,cr2,cr3,dr2,dr3,ti2,tr2;
    t0=l1*ido;

    t1=0;
    t2=t0<<1;
    t3=ido<<1;
    t4=ido+(ido<<1);
    t5=0;
    for(k=0;k<l1;k++)
    {
        tr2=cc[t3-1]+cc[t3-1];
        cr2=cc[t5]+(taur*tr2);
        ch[t1]=cc[t5]+tr2;
        ci3=taui*(cc[t3]+cc[t3]);
        ch[t1+t0]=cr2-ci3;
        ch[t1+t2]=cr2+ci3;
        t1+=ido;
        t3+=t4;
        t5+=t4;
    }

    if(ido==1)return;

    t1=0;
    t3=ido<<1;
    for(k=0;k<l1;k++)
    {
        t7=t1+(t1<<1);
        t6=(t5=t7+t3);
        t8=t1;
        t10=(t9=t1+t0)+t0;

        for(i=2;i<ido;i+=2)
        {
            t5+=2;
            t6-=2;
            t7+=2;
            t8+=2;
            t9+=2;
            t10+=2;
            tr2=cc[t5-1]+cc[t6-1];
            cr2=cc[t7-1]+(taur*tr2);
            ch[t8-1]=cc[t7-1]+tr2;
            ti2=cc[t5]-cc[t6];
            ci2=cc[t7]+(taur*ti2);
            ch[t8]=cc[t7]+ti2;
            cr3=taui*(cc[t5-1]-cc[t6-1]);
            ci3=taui*(cc[t5]+cc[t6]);
            dr2=cr2-ci3;
            dr3=cr2+ci3;
            di2=ci2+cr3;
            di3=ci2-cr3;
            ch[t9-1]=wa1[i-2]*dr2-wa1[i-1]*di2;
            ch[t9]=wa1[i-2]*di2+wa1[i-1]*dr2;
            ch[t10-1]=wa2[i-2]*dr3-wa2[i-1]*di3;
            ch[t10]=wa2[i-2]*di3+wa2[i-1]*dr3;
        }
        t1+=ido;
    }
}

static void dradb4(OsInt32 ido,OsInt32 l1,OsFloat *cc,OsFloat *ch,OsFloat *wa1,OsFloat *wa2,OsFloat *wa3)
{
    static OsFloat sqrt2=1.414213562373095f;
    OsInt32 i,k,t0,t1,t2,t3,t4,t5,t6,t7,t8;
    OsFloat ci2,ci3,ci4,cr2,cr3,cr4,ti1,ti2,ti3,ti4,tr1,tr2,tr3,tr4;
    t0=l1*ido;

    t1=0;
    t2=ido<<2;
    t3=0;
    t6=ido<<1;
    for(k = 0; k < l1; k++)
    {
        t4=t3+t6;
        t5=t1;
        tr3=cc[t4-1]+cc[t4-1];
        tr4=cc[t4]+cc[t4]; 
        tr1=cc[t3]-cc[(t4+=t6)-1];
        tr2=cc[t3]+cc[t4-1];
        ch[t5]=tr2+tr3;
        ch[t5+=t0]=tr1-tr4;
        ch[t5+=t0]=tr2-tr3;
        ch[t5+=t0]=tr1+tr4;
        t1+=ido;
        t3+=t2;
    }
    
    if(ido < 2) return;
    if(ido == 2) goto L105;

    t1=0;
    for(k=0;k<l1;k++)
    {
        t5=(t4=(t3=(t2=t1<<2)+t6))+t6;
        t7=t1;
        for(i=2;i<ido;i+=2)
        {
            t2+=2;
            t3+=2;
            t4-=2;
            t5-=2;
            t7+=2;
            ti1=cc[t2]+cc[t5];
            ti2=cc[t2]-cc[t5];
            ti3=cc[t3]-cc[t4];
            tr4=cc[t3]+cc[t4];
            tr1=cc[t2-1]-cc[t5-1];
            tr2=cc[t2-1]+cc[t5-1];
            ti4=cc[t3-1]-cc[t4-1];
            tr3=cc[t3-1]+cc[t4-1];
            ch[t7-1]=tr2+tr3;
            cr3=tr2-tr3;
            ch[t7]=ti2+ti3;
            ci3=ti2-ti3;
            cr2=tr1-tr4;
            cr4=tr1+tr4;
            ci2=ti1+ti4;
            ci4=ti1-ti4;

            ch[(t8=t7+t0)-1]=wa1[i-2]*cr2-wa1[i-1]*ci2;
            ch[t8]=wa1[i-2]*ci2+wa1[i-1]*cr2;
            ch[(t8+=t0)-1]=wa2[i-2]*cr3-wa2[i-1]*ci3;
            ch[t8]=wa2[i-2]*ci3+wa2[i-1]*cr3;
            ch[(t8+=t0)-1]=wa3[i-2]*cr4-wa3[i-1]*ci4;
            ch[t8]=wa3[i-2]*ci4+wa3[i-1]*cr4;
        }
        t1+=ido;
    }
    if(ido%2 == 1)return;
L105:
    t1=ido;
    t2=ido<<2;
    t3=ido-1;
    t4=ido+(ido<<1);
    for(k=0;k<l1;k++)
    {
        t5=t3;
        ti1=cc[t1]+cc[t4];
        ti2=cc[t4]-cc[t1];
        tr1=cc[t1-1]-cc[t4-1];
        tr2=cc[t1-1]+cc[t4-1];
        ch[t5]=tr2+tr2;
        ch[t5+=t0]=sqrt2*(tr1-ti1);
        ch[t5+=t0]=ti2+ti2;
        ch[t5+=t0]=-sqrt2*(tr1+ti1);

        t3+=ido;
        t1+=t2;
        t4+=t2;
    }
}

static void dradbg(OsInt32 ido,OsInt32 ip,OsInt32 l1,OsInt32 idl1,OsFloat *cc,OsFloat *c1,OsFloat *c2,OsFloat *ch,OsFloat *ch2,OsFloat *wa)
{
    static OsFloat tpi = 6.283185307179586f;
    OsInt32 idij,ipph,i,j,k,l,ik,is,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12;
    OsFloat dc2,ai1,ai2,ar1,ar2,ds2;
    OsInt32 nbd;
    OsFloat dcp,arg,dsp,ar1h,ar2h;
    OsInt32 ipp2;

    t10 = ip*ido;
    t0 = l1*ido;
    arg = tpi/(OsFloat)ip;
    dcp = cosf(arg);
    dsp = sinf(arg);
    nbd = (ido-1)>>1;
    ipp2 = ip;
    ipph = (ip+1)>>1;
    if(ido < l1) goto L103;

    t1 = 0;
    t2 = 0;
    for(k = 0; k < l1; k++)
    {
        t3=t1;
        t4=t2;
        for(i=0;i<ido;i++)
        {
            ch[t3]=cc[t4];
            t3++;
            t4++;
        }
        t1+=ido;
        t2+=t10;
    }
    goto L106;

L103:
    t1=0;
    for(i=0;i<ido;i++)
    {
        t2=t1;
        t3=t1;
        for(k=0;k<l1;k++)
        {
            ch[t2]=cc[t3];
            t2+=ido;
            t3+=t10;
        }
        t1++;
    }

L106:
    t1=0;
    t2=ipp2*t0;
    t7=(t5=ido<<1);
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;
        t6=t5;
        for(k=0;k<l1;k++)
        {
            ch[t3]=cc[t6-1]+cc[t6-1];
            ch[t4]=cc[t6]+cc[t6];
            t3+=ido;
            t4+=ido;
            t6+=t10;
        }
        t5+=t7;
    }

    if (ido == 1)goto L116;
    if(nbd<l1)goto L112;

    t1=0;
    t2=ipp2*t0;
    t7=0;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;

        t7+=(ido<<1);
        t8=t7;
        for(k=0;k<l1;k++)
        {
            t5=t3;
            t6=t4;
            t9=t8;
            t11=t8;
            for(i=2;i<ido;i+=2)
            {
                t5+=2;
                t6+=2;
                t9+=2;
                t11-=2;
                ch[t5-1]=cc[t9-1]+cc[t11-1];
                ch[t6-1]=cc[t9-1]-cc[t11-1];
                ch[t5]=cc[t9]-cc[t11];
                ch[t6]=cc[t9]+cc[t11];
            }
            t3+=ido;
            t4+=ido;
            t8+=t10;
        }
    }
    goto L116;

L112:
    t1=0;
    t2=ipp2*t0;
    t7=0;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;
        t7+=(ido<<1);
        t8=t7;
        t9=t7;
        for(i=2;i<ido;i+=2)
        {
            t3+=2;
            t4+=2;
            t8+=2;
            t9-=2;
            t5=t3;
            t6=t4;
            t11=t8;
            t12=t9;
            for(k=0;k<l1;k++)
            {
                ch[t5-1]=cc[t11-1]+cc[t12-1];
                ch[t6-1]=cc[t11-1]-cc[t12-1];
                ch[t5]=cc[t11]-cc[t12];
                ch[t6]=cc[t11]+cc[t12];
                t5+=ido;
                t6+=ido;
                t11+=t10;
                t12+=t10;
            }
        }
    }

L116:
    ar1=1.f;
    ai1=0.f;
    t1=0;
    t9=(t2=ipp2*idl1);
    t3=(ip-1)*idl1;
    for(l=1;l<ipph;l++)
    {
        t1+=idl1;
        t2-=idl1;

        ar1h=dcp*ar1-dsp*ai1;
        ai1=dcp*ai1+dsp*ar1;
        ar1=ar1h;
        t4=t1;
        t5=t2;
        t6=0;
        t7=idl1;
        t8=t3;
        for(ik=0;ik<idl1;ik++)
        {
            c2[t4++]=ch2[t6++]+ar1*ch2[t7++];
            c2[t5++]=ai1*ch2[t8++];
        }
        dc2=ar1;
        ds2=ai1;
        ar2=ar1;
        ai2=ai1;

        t6=idl1;
        t7=t9-idl1;
        for(j=2;j<ipph;j++)
        {
            t6+=idl1;
            t7-=idl1;
            ar2h=dc2*ar2-ds2*ai2;
            ai2=dc2*ai2+ds2*ar2;
            ar2=ar2h;
            t4=t1;
            t5=t2;
            t11=t6;
            t12=t7;
            for(ik=0;ik<idl1;ik++)
            {
                c2[t4++]+=ar2*ch2[t11++];
                c2[t5++]+=ai2*ch2[t12++];
            }
        }
    }

    t1=0;
    for(j=1;j<ipph;j++)
    {
        t1+=idl1;
        t2=t1;
        for(ik=0;ik<idl1;ik++)ch2[ik]+=ch2[t2++];
    }

    t1=0;
    t2=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;
        for(k=0;k<l1;k++)
        {
            ch[t3]=c1[t3]-c1[t4];
            ch[t4]=c1[t3]+c1[t4];
            t3+=ido;
            t4+=ido;
        }
    }

    if(ido==1)goto L132;
    if(nbd<l1)goto L128;

    t1=0;
    t2=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;
        for(k=0;k<l1;k++)
        {
            t5=t3;
            t6=t4;
            for(i=2;i<ido;i+=2)
            {
                t5+=2;
                t6+=2;
                ch[t5-1]=c1[t5-1]-c1[t6];
                ch[t6-1]=c1[t5-1]+c1[t6];
                ch[t5]=c1[t5]+c1[t6-1];
                ch[t6]=c1[t5]-c1[t6-1];
            }
            t3+=ido;
            t4+=ido;
        }
    }
    goto L132;

L128:
    t1=0;
    t2=ipp2*t0;
    for(j=1;j<ipph;j++)
    {
        t1+=t0;
        t2-=t0;
        t3=t1;
        t4=t2;
        for(i=2;i<ido;i+=2)
        {
            t3+=2;
            t4+=2;
            t5=t3;
            t6=t4;
            for(k=0;k<l1;k++)
            {
                ch[t5-1]=c1[t5-1]-c1[t6];
                ch[t6-1]=c1[t5-1]+c1[t6];
                ch[t5]=c1[t5]+c1[t6-1];
                ch[t6]=c1[t5]-c1[t6-1];
                t5+=ido;
                t6+=ido;
            }
        }
    }

L132:
    if(ido==1)return;

    for(ik=0;ik<idl1;ik++)c2[ik]=ch2[ik];

    t1=0;
    for(j=1;j<ip;j++)
    {
        t2=(t1+=t0);
        for(k=0;k<l1;k++)
        {
            c1[t2]=ch[t2];
            t2+=ido;
        }
    }

    if(nbd>l1)goto L139;

    is= -ido-1;
    t1=0;
    for(j=1;j<ip;j++)
    {
        is+=ido;
        t1+=t0;
        idij=is;
        t2=t1;
        for(i=2;i<ido;i+=2)
        {
            t2+=2;
            idij+=2;
            t3=t2;
            for(k=0;k<l1;k++)
            {
                c1[t3-1]=wa[idij-1]*ch[t3-1]-wa[idij]*ch[t3];
                c1[t3]=wa[idij-1]*ch[t3]+wa[idij]*ch[t3-1];
                t3+=ido;
            }
        }
    }
    return;

L139:
    is= -ido-1;
    t1=0;
    for(j=1;j<ip;j++)
    {
        is+=ido;
        t1+=t0;
        t2=t1;
        for(k=0;k<l1;k++)
        {
            idij=is;
            t3=t2;
            for(i=2;i<ido;i+=2)
            {
                idij+=2;
                t3+=2;
                c1[t3-1]=wa[idij-1]*ch[t3-1]-wa[idij]*ch[t3];
                c1[t3]=wa[idij-1]*ch[t3]+wa[idij]*ch[t3-1];
            }
            t2+=ido;
        }
    }
}

static void drftb1(OsInt32 n,OsFloat *c,OsFloat *ch,OsFloat *wa,OsInt32 *ifac)
{
    OsInt32 i,k1,l1,l2;
    OsInt32 na;
    OsInt32 nf,ip,iw,ix2,ix3,ido,idl1;

    nf  = ifac[1];
    na  = 0;
    l1  = 1;
    iw  = 1;

    for(k1=0;k1<nf;k1++)
    {
        ip = ifac[k1 + 2];
        l2 = ip*l1;
        ido = n/l2;
        idl1= ido*l1;
        if(ip != 4) goto L103;
        ix2 = iw+ido;
        ix3 = ix2+ido;

        if(na != 0)
            dradb4(ido,l1,ch,c,wa+iw-1,wa+ix2-1,wa+ix3-1);
        else
            dradb4(ido,l1,c,ch,wa+iw-1,wa+ix2-1,wa+ix3-1);
        na = 1-na;
        goto L115;
L103:
        if(ip!=2) goto L106;
        if(na!=0)
            dradb2(ido,l1,ch,c,wa+iw-1);
        else
            dradb2(ido,l1,c,ch,wa+iw-1);
        na = 1-na;
        goto L115;
L106:
        if(ip != 3) goto L109;
        ix2 = iw+ido;
        if(na != 0)
            dradb3(ido,l1,ch,c,wa+iw-1,wa+ix2-1);
        else
            dradb3(ido,l1,c,ch,wa+iw-1,wa+ix2-1);
        na = 1-na;
        goto L115;
L109:
        if(na != 0)
            dradbg(ido,ip,l1,idl1,ch,ch,ch,c,c,wa+iw-1);
        else
            dradbg(ido,ip,l1,idl1,c,c,c,ch,ch,wa+iw-1);
        if(ido == 1) na = 1-na;
L115:
        l1 = l2;
        iw += (ip-1)*ido;
    }
    if(na == 0) return;
    for(i = 0; i < n; i++) c[i] = ch[i];
}

void kiss_drft_forward(kft_t *l,OsFloat *data)
{
    if(l->n == 1) return;
    drftf1(l->n,data,l->trigcache,l->trigcache+l->n,l->splitcache);
}

void kiss_drft_backward(kft_t *l,OsFloat *data)
{
    if (l->n == 1) return;
    drftb1(l->n,data,l->trigcache,l->trigcache+l->n,l->splitcache);
}

void kiss_drft_init(kft_t *l,OsInt32 n)
{
    l->n = n;
    l->trigcache = (OsFloat*)icm_alloc(3*n*sizeof(*l->trigcache),1);
    l->splitcache = (OsInt32*)icm_alloc(32*sizeof(*l->splitcache),1);
    fdrffti(n,l->trigcache,l->splitcache);
}

void kiss_drft_clear(kft_t *l)
{
    if(l)
    {
        if(l->trigcache)
            icm_free(l->trigcache);
        if(l->splitcache)
            icm_free(l->splitcache);
    }
}

/***************************************** cdft ******************************************/

struct cdft_t
{
    OsInt32 bits;
    OsFloat *sin_table;
    OsFloat *cos_table;
};

cdft_t* icm_cdft_init(OsInt32 M)
{
    OsInt32 bits = ceil(log(1.0*M)/log(2.0));
    OsFloat twopi_n,twopi_num;
    OsInt32 i;
    cdft_t *fft;

    if(M != pow(2.0,bits)) return 0;
    twopi_n     = (atan(1.0)*8.0) / M;
    twopi_num   = (atan(1.0)*8.0) / (M-1);

    fft = (cdft_t*)icm_alloc(sizeof(cdft_t),1);
    if (!fft) return NULL;

    fft->bits       = bits;
    fft->sin_table  = (OsFloat*)icm_alloc(M,sizeof(OsFloat));
    fft->cos_table  = (OsFloat*)icm_alloc(M,sizeof(OsFloat));

    for (i = 0; i < M; i++)
    {
        fft->sin_table[i] = sin(i*twopi_n);
        fft->cos_table[i] = cos(i*twopi_n);
    }
    return fft;
}

void icm_cdft_uninit(cdft_t *fft)
{
    if (fft)
    {
        icm_free(fft->sin_table);
        icm_free(fft->cos_table);
        icm_free(fft);
    }
}

static const OsUInt32 reverse(OsUInt32 val,OsInt32 bits)
{
    OsUInt32 retn = 0;
    while (bits--)
    {
        retn <<= 1;
        retn |= (val & 1);
        val >>= 1;
    }
    return retn;
}

static OsFloat fft_amp(OsInt32 n,complex x[],OsInt32 bits)
{
    n = reverse(n,bits);
    return (hypot(x[n].re,x[n].im));
}

static void fft_scale(complex x[],OsInt32 bits)
{
    OsInt32 i;
    OsInt32 samples = (1 << bits);

    for (i = 0; i < samples; i++)
    {
        x[i].re /= samples;
        x[i].im /= samples;
    }
}

void icm_cdft_fft(cdft_t *fft,complex x[])
{
    register OsInt32 loop,loop1,loop2;
    OsUInt32    i1;
    OsInt32     i2,i3,i4,y;
    OsFloat     a1,a2,b1,b2,z1,z2;
    OsInt32     bits = fft->bits;
    OsInt32     samples;

    complex   t;

    i1 = (1 << bits) / 2;
    i2 = 1;

    for (loop = 0; loop < bits; loop++)
    {
        i3 = 0;
        i4 = i1;

        for (loop1 = 0; loop1 < i2; loop1++)
        {
            y  =  reverse(i3/(OsInt32)i1,bits);
            z1 =  fft->cos_table[y];
            z2 = -fft->sin_table[y];

            for (loop2 = i3; loop2 < i4; loop2++)
            {
                a1 = x[loop2].re;
                a2 = x[loop2].im;

                b1 = z1 * x[loop2+i1].re - z2 * x[loop2+i1].im;
                b2 = z2 * x[loop2+i1].re + z1 * x[loop2+i1].im;

                x[loop2].re = a1 + b1;
                x[loop2].im = a2 + b2;

                x[loop2+i1].re = a1 - b1;
                x[loop2+i1].im = a2 - b2;
            }
            i3 += (i1 << 1);
            i4 += (i1 << 1);
        }
        i1 >>= 1;
        i2 <<= 1;
    }
    samples = (1 << bits)/2;
    for(loop = 0; loop < samples; loop++)
    {
        OsUInt32 r = reverse(loop,bits);
        if(r == loop) continue;
        t = x[loop];
        x[loop] = x[r];
        x[r] = t;
    }
}

void icm_cdft_ifft(cdft_t *fft,complex x[])
{
    register OsInt32 loop,loop1,loop2;
    OsUInt32    i1;
    OsInt32     i2,i3,i4,y;
    OsFloat     a1,a2,b1,b2,z1,z2;
    OsInt32     bits = fft->bits;
    OsInt32     samples;
    complex     t;

    i1 = (1 << bits) / 2;
    i2 = 1;

    for (loop = 0; loop < bits; loop++)
    {
        i3 = 0;
        i4 = i1;

        for (loop1 = 0; loop1 < i2; loop1++)
        {
            y  =  reverse(i3/(OsInt32)i1,bits);
            z1 =  fft->cos_table[y];
            z2 =  fft->sin_table[y];

            for (loop2 = i3; loop2 < i4; loop2++)
            {
                a1 = x[loop2].re;
                a2 = x[loop2].im;

                b1 = z1 * x[loop2+i1].re - z2 * x[loop2+i1].im;
                b2 = z2 * x[loop2+i1].re + z1 * x[loop2+i1].im;

                x[loop2].re = a1 + b1;
                x[loop2].im = a2 + b2;

                x[loop2+i1].re = a1 - b1;
                x[loop2+i1].im = a2 - b2;
            }
            i3 += (i1 << 1);
            i4 += (i1 << 1);
        }
        i1 >>= 1;
        i2 <<= 1;
    }
    samples = (1 << bits)/2;
    for(loop = 0; loop < samples; loop++)
    {
        OsUInt32 r = reverse(loop,bits);
        if(r == loop) continue;
        t = x[loop];
        x[loop] = x[r];
        x[r] = t;
    }
    fft_scale(x,bits);
}

/*
//    OsFloat x[N] = {1.0f,0.7071f,0.0f,-0.7071f,-1.0f,-0.7071f,-0.0f,0.7071f};
    OsFloat x[N] = {1.0f,-0.9239f,0.7071f,-0.3827f,-0.0000f,0.3827f,-0.7071f,0.9239f};
    OsInt32 i;

    // rdft
    fft_t *fft = icm_fft_init(N);
    complex X[N/2+1];
    icm_fft(fft,x,X);
    printf("\n");
    for(i = 0; i < N/2+1; i++)
    {
        printf("re=%f,im=%f\n",X[i].re,X[i].im);
    }
    printf("\n");
    OsFloat x2[N] = {0};
    icm_ifft(fft,X,x2);
    for(i = 0; i < N; i++)
    {
        printf("%f\n",x2[i]);
    }
    icm_fft_uninit(fft);

    fft = icm_fft_init2(N);
    icm_fft2(fft,x);
    printf("\nre=%f,im=%f\n",x[0],0.0);
    for(i = 1; i < N/2; i++)
    {
        printf("re=%f,im=%f\n",x[2*i],x[2*i+1]);
    }
    printf("re=%f,im=%f\n",x[1],0.0);
    printf("\n");
    icm_ifft2(fft,x);
    for(i = 0; i < N; i++)
    {
        printf("%f\n",x[i]);
    }
    printf("\n");
    icm_fft_uninit2(fft);

    // speex fft
    void *table = kiss_fft_init(N);
    OsFloat x3[N] = {0};
    kiss_fft(table,x,x3);
    printf("\n");
    printf("re=%f,im=%f\n",x3[0],0.0);
    for(i = 1; i < N/2; i++)
    {
        printf("re=%f,im=%f\n",x3[2*i-1],x3[2*i]);
    }
    printf("re=%f,im=%f\n",x3[N-1],0.0);
    printf("\n");

    OsFloat x4[N] = {0};
    kiss_ifft(fft,x3,x4);
    for(i = 0; i < N; i++)
    {
        printf("%f\n",x4[i]);
    }
    printf("\n");
    kiss_fft_uninit(table);
*/