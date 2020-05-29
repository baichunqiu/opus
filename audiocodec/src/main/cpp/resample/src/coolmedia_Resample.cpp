
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "coolmedia_Resample.h"

//////////////////////////////////////////////////////////////////////
/*
 * Taken from stddefs.h
 */
#ifndef PI
#define PI (3.14159265358979232846)
#endif

#ifndef PI2
#define PI2 (6.28318530717958465692)
#endif

#define D2R (0.01745329348)          /* (2*pi)/360 */
#define R2D (57.29577951)            /* 360/(2*pi) */

#ifndef MAX
#define MAX(x,y) ((x)>(y) ?(x):(y))
#endif
#ifndef MIN
#define MIN(x,y) ((x)<(y) ?(x):(y))
#endif

#ifndef ABS
#define ABS(x) ((x)<0 ?(-(x)):(x))
#endif

#ifndef SGN
#define SGN(x) ((x)<0 ? (-1) : ((x)==0 ? (0) : (1)))
#endif

typedef char           RES_BOOL;
typedef short          RES_HWORD;
typedef int            RES_WORD;
typedef unsigned short RES_UHWORD;
typedef unsigned int   RES_UWORD;

#define MAX_HWORD (32767)
#define MIN_HWORD (-32768)

#ifdef DEBUG
#define INLINE
#else
#define INLINE inline
#endif

/*
 * Taken from resample.h
 *
 * The configuration constants below govern
 * the number of bits in the input sample and filter coefficients, the
 * number of bits to the right of the binary-point for fixed-point math, etc.
 *
 */

/* Conversion constants */
#define Nhc     8
#define Na      7
#define Np      (Nhc+Na)
#define Npc     (1<<Nhc)
#define Amask   ((1<<Na)-1)
#define Pmask   ((1<<Np)-1)
#define Nh      16
#define Nb      16
#define Nhxn    14
#define Nhg     (Nh-Nhxn)
#define NLpScl  13

/* Description of constants:
 *
 * Npc - is the number of look-up values available for the lowpass filter
 *    between the beginning of its impulse response and the "cutoff time"
 *    of the filter.  The cutoff time is defined as the reciprocal of the
 *    lowpass-filter cut off frequence in Hz.  For example, if the
 *    lowpass filter were a sinc function, Npc would be the index of the
 *    impulse-response lookup-table corresponding to the first zero-
 *    crossing of the sinc function.  (The inverse first zero-crossing
 *    time of a sinc function equals its nominal cutoff frequency in Hz.)
 *    Npc must be a power of 2 due to the details of the current
 *    implementation. The default value of 512 is sufficiently high that
 *    using linear interpolation to fill in between the table entries
 *    gives approximately 16-bit accuracy in filter coefficients.
 *
 * Nhc - is log base 2 of Npc.
 *
 * Na - is the number of bits devoted to linear interpolation of the
 *    filter coefficients.
 *
 * Np - is Na + Nhc, the number of bits to the right of the binary point
 *    in the integer "time" variable. To the left of the point, it indexes
 *    the input array (X), and to the right, it is interpreted as a number
 *    between 0 and 1 sample of the input X.  Np must be less than 16 in
 *    this implementation.
 *
 * Nh - is the number of bits in the filter coefficients. The sum of Nh and
 *    the number of bits in the input data (typically 16) cannot exceed 32.
 *    Thus Nh should be 16.  The largest filter coefficient should nearly
 *    fill 16 bits (32767).
 *
 * Nb - is the number of bits in the input data. The sum of Nb and Nh cannot
 *    exceed 32.
 *
 * Nhxn - is the number of bits to right shift after multiplying each input
 *    sample times a filter coefficient. It can be as great as Nh and as
 *    small as 0. Nhxn = Nh-2 gives 2 guard bits in the multiply-add
 *    accumulation.  If Nhxn=0, the accumulation will soon overflow 32 bits.
 *
 * Nhg - is the number of guard bits in mpy-add accumulation (equal to Nh-Nhxn)
 *
 * NLpScl - is the number of bits allocated to the unity-gain normalization
 *    factor.  The output of the lowpass filter is multiplied by LpScl and
 *    then right-shifted NLpScl bits. To avoid overflow, we must have
 *    Nb+Nhg+NLpScl < 32.
 */


#ifdef _MSC_VER
#pragma warning(push, 3)
//#pragma warning(disable: 4245)  // Conversion from uint to ushort
#pragma warning(disable: 4244)  // Conversion from double to uint
#pragma warning(disable: 4146)  // unary minus operator applied to unsigned type, result still unsigned
#pragma warning(disable: 4761)  // integral size mismatch in argument; conversion supplied
#endif

#if 1
#include "smallfilter.h"
#else
#define SMALL_FILTER_NMULT  0
#define SMALL_FILTER_SCALE  0
#define SMALL_FILTER_NWING  0
#define SMALL_FILTER_IMP    NULL
#define SMALL_FILTER_IMPD   NULL
#endif

#if 1
#include "largefilter.h"
#else
#define LARGE_FILTER_NMULT  0
#define LARGE_FILTER_SCALE  0
#define LARGE_FILTER_NWING  0
#define LARGE_FILTER_IMP    NULL
#define LARGE_FILTER_IMPD   NULL
#endif

#undef  INLINE
#define INLINE
#define HAVE_FILTER 0

#ifndef NULL
#define NULL    0
#endif

static INLINE RES_HWORD WordToHword(RES_WORD v, int scl)
{
    RES_HWORD out;
    RES_WORD llsb = (1<<(scl-1));
    v += llsb;          /* round */
    v >>= scl;
    if(v>MAX_HWORD)
    {
        v = MAX_HWORD;
    }
    else if (v < MIN_HWORD)
    {
        v = MIN_HWORD;
    }
    out = (RES_HWORD)v;
    return out;
}

// Sampling rate conversion using linear interpolation for maximum speed.
static int SrcLinear(const RES_HWORD X[],RES_HWORD Y[],double pFactor,RES_UHWORD nx)
{
    RES_HWORD iconst;
    RES_UWORD time = 0;
    const RES_HWORD *xp;
    RES_HWORD *Ystart, *Yend;
    RES_WORD v,x1,x2;

    double dt;                              /* Step through input signal */
    RES_UWORD dtb;                          /* Fixed-point version of Dt */
    RES_UWORD endTime;                      /* When time reaches EndTime, return to user */

    dt = 1.0/pFactor;                       /* Output sampling period */
    dtb = dt*(1<<Np) + 0.5;                 /* Fixed-point representation */

    Ystart = Y;
    Yend = Ystart + (unsigned)(nx * pFactor);
    endTime = time + (1<<Np)*(RES_WORD)nx;
    while (time < endTime)
    {
        iconst = (time) & Pmask;
        xp = &X[(time)>>Np];                /* Ptr to current input sample */
        x1 = *xp++;
        x2 = *xp;
        x1 *= ((1<<Np)-iconst);
        x2 *= iconst;
        v = x1 + x2;
        *Y++ = WordToHword(v,Np);           /* Deposit output */
        time += dtb;                        /* Move to next sample by time increment */
    }
    return (Y - Ystart);                    /* Return number of output samples */
}

static RES_WORD FilterUp(const RES_HWORD Imp[],const RES_HWORD ImpD[],
                        RES_UHWORD Nwing,RES_BOOL Interp,
                        const RES_HWORD *Xp,RES_HWORD Ph,RES_HWORD Inc)
{
    const RES_HWORD *Hp;
    const RES_HWORD *Hdp = NULL;
    const RES_HWORD *End;
    RES_HWORD a = 0;
    RES_WORD v, t;

    v = 0;
    Hp = &Imp[Ph>>Na];
    End = &Imp[Nwing];
    if(Interp)
    {
        Hdp = &ImpD[Ph>>Na];
        a = Ph & Amask;
    }
    if(Inc == 1)                            /* If doing right wing...              */
    {                                       /* ...drop extra coeff, so when Ph is  */
        End--;                              /*    0.5, we don't do too many mult's */
        if (Ph == 0)                        /* If the phase is zero...           */
        {                                   /* ...then we've already skipped the */
            Hp += Npc;                      /*    first sample, so we must also  */
            Hdp += Npc;                     /*    skip ahead in Imp[] and ImpD[] */
        }
    }
    if(Interp)
    {
        while (Hp < End)
        {
            t = *Hp;                /* Get filter coeff */
            t += (((RES_WORD)*Hdp)*a)>>Na; /* t is now interp'd filter coeff */
            Hdp += Npc;             /* Filter coeff differences step */
            t *= *Xp;               /* Mult coeff by input sample */
            if (t & (1<<(Nhxn-1)))  /* Round, if needed */
                t += (1<<(Nhxn-1));
            t >>= Nhxn;             /* Leave some guard bits, but come back some */
            v += t;                 /* The filter output */
            Hp += Npc;              /* Filter coeff step */
            Xp += Inc;              /* Input signal step. NO CHECK ON BOUNDS */
        }
    }
    else
    {
        while (Hp < End)
        {
            t = *Hp;                /* Get filter coeff */
            t *= *Xp;               /* Mult coeff by input sample */
            if (t & (1<<(Nhxn-1)))  /* Round, if needed */
                t += (1<<(Nhxn-1));
            t >>= Nhxn;             /* Leave some guard bits, but come back some */
            v += t;	                /* The filter output */
            Hp += Npc;              /* Filter coeff step */
            Xp += Inc;              /* Input signal step. NO CHECK ON BOUNDS */
        }
    }
    return(v);
}


static RES_WORD FilterUD(const RES_HWORD Imp[],const RES_HWORD ImpD[],
                        RES_UHWORD Nwing,RES_BOOL Interp,
                        const RES_HWORD *Xp,RES_HWORD Ph,RES_HWORD Inc,RES_UHWORD dhb)
{
    RES_HWORD a;
    const RES_HWORD *Hp, *Hdp, *End;
    RES_WORD v, t;
    RES_UWORD Ho;

    v = 0;
    Ho = (Ph*(RES_UWORD)dhb)>>Np;
    End = &Imp[Nwing];
    if (Inc == 1)                           /* If doing right wing...              */
    {                                       /* ...drop extra coeff, so when Ph is  */
        End--;                          /*    0.5, we don't do too many mult's */
        if (Ph == 0)                    /* If the phase is zero...           */
            Ho += dhb;              /* ...then we've already skipped the */
    }                                       /*    first sample, so we must also  */
                                                /*    skip ahead in Imp[] and ImpD[] */
    if(Interp)
    {
        while ((Hp = &Imp[Ho>>Na]) < End)
        {
            t = *Hp;                /* Get IR sample */
            Hdp = &ImpD[Ho>>Na];    /* get interp (lower Na) bits from diff table*/
            a = Ho & Amask;         /* a is logically between 0 and 1 */
            t += (((RES_WORD)*Hdp)*a)>>Na;  /* t is now interp'd filter coeff */
            t *= *Xp;               /* Mult coeff by input sample */
            if (t & 1<<(Nhxn-1))    /* Round, if needed */
                t += 1<<(Nhxn-1);
            t >>= Nhxn;             /* Leave some guard bits, but come back some */
            v += t;                 /* The filter output */
            Ho += dhb;              /* IR step */
            Xp += Inc;              /* Input signal step. NO CHECK ON BOUNDS */
        }
    }
    else
    {
        while ((Hp = &Imp[Ho>>Na]) < End)
        {
            t = *Hp;                /* Get IR sample */
            t *= *Xp;               /* Mult coeff by input sample */
                        if (t & 1<<(Nhxn-1))    /* Round, if needed */
                                t += 1<<(Nhxn-1);
                        t >>= Nhxn;             /* Leave some guard bits, but come back some */
                        v += t;                 /* The filter output */
                        Ho += dhb;              /* IR step */
                        Xp += Inc;              /* Input signal step. NO CHECK ON BOUNDS */
                }
        }
        return(v);
}

/* Sampling rate up-conversion only subroutine;
 * Slightly faster than down-conversion;
 */
static int SrcUp(const RES_HWORD X[], RES_HWORD Y[], double pFactor,
                RES_UHWORD nx, RES_UHWORD pNwing, RES_UHWORD pLpScl,
                const RES_HWORD pImp[], const RES_HWORD pImpD[], RES_BOOL Interp)
{
    const RES_HWORD *xp;
    RES_HWORD *Ystart, *Yend;
    RES_WORD v;

    double dt;                          /* Step through input signal */
    RES_UWORD dtb;                      /* Fixed-point version of Dt */
    RES_UWORD time = 0;
    RES_UWORD endTime;                  /* When time reaches EndTime, return to user */

    dt = 1.0/pFactor;                   /* Output sampling period */
    dtb = dt*(1<<Np) + 0.5;             /* Fixed-point representation */

    Ystart = Y;
    Yend = Ystart + (unsigned)(nx * pFactor);
    endTime = time + (1<<Np)*(RES_WORD)nx;
    while(time < endTime)
    {
        xp = &X[time>>Np];              /* Ptr to current input sample */
        /* Perform left-wing inner product */
        v = 0;
        v = FilterUp(pImp, pImpD, pNwing, Interp, xp, (RES_HWORD)(time&Pmask),-1);
        /* Perform right-wing inner product */
        v += FilterUp(pImp, pImpD, pNwing, Interp, xp+1,  (RES_HWORD)((-time)&Pmask),1);
        v >>= Nhg;                      /* Make guard bits */
        v *= pLpScl;                    /* Normalize for unity filter gain */
        *Y++ = WordToHword(v,NLpScl);   /* strip guard bits, deposit output */
        time += dtb;                    /* Move to next sample by time increment */
    }
    return (Y - Ystart);                /* Return the number of output samples */
}

/* Sampling rate conversion subroutine */
static int SrcUD(const RES_HWORD X[],RES_HWORD Y[],double pFactor,
                RES_UHWORD nx,RES_UHWORD pNwing,RES_UHWORD pLpScl,
                const RES_HWORD pImp[],const RES_HWORD pImpD[],RES_BOOL Interp)
{
    const RES_HWORD *xp;
    RES_HWORD *Ystart, *Yend;
    RES_WORD v;

    double dh;                          /* Step through filter impulse response */
    double dt;                          /* Step through input signal */
    RES_UWORD time = 0;
    RES_UWORD endTime;                  /* When time reaches EndTime, return to user */
    RES_UWORD dhb, dtb;                 /* Fixed-point versions of Dh,Dt */

    dt = 1.0/pFactor;                   /* Output sampling period */
    dtb = dt*(1<<Np) + 0.5;             /* Fixed-point representation */

    dh = MIN(Npc, pFactor*Npc);         /* Filter sampling period */
    dhb = dh*(1<<Na) + 0.5;             /* Fixed-point representation */

    Ystart = Y;
    Yend = Ystart + (unsigned)(nx * pFactor);
    endTime = time + (1<<Np)*(RES_WORD)nx;
    while(time < endTime)
    {
        xp = &X[time>>Np];              /* Ptr to current input sample */
        /* Perform left-wing inner product */
        v = FilterUD(pImp,pImpD,pNwing,Interp,xp,(RES_HWORD)(time&Pmask),-1,dhb);
        /* Perform right-wing inner product */
        v += FilterUD(pImp,pImpD,pNwing,Interp,xp+1,(RES_HWORD)((-time)&Pmask),1,dhb);
        v >>= Nhg;                      /* Make guard bits */
        v *= pLpScl;                    /* Normalize for unity filter gain */
        *Y++ = WordToHword(v,NLpScl);   /* strip guard bits, deposit output */
        time += dtb;                    /* Move to next sample by time increment */
    }
    return (Y - Ystart);                /* Return the number of output samples */
}

Resample::Resample()
    : m_nInSampleRate(0)
    , m_nOutSampleRate(0)
    , m_nSamplesPerFrame(0)
{
    m_factor        = 1.0;              /* Conversion factor = nOutSampleRate / nInSampleRate.  */
    m_large_filter  = false;            /* Large filter?			    */
    m_high_quality  = true;             /* Not fast?				    */
    m_xoff          = 0;                /* History and lookahead size, in samples   */
    m_buffer        = NULL;             /* Input buffer.			    */
}

Resample::~Resample()
{
    Uninit();
}

int Resample::Init(int nInSampleRate,int nOutSampleRate,int nSamplesPerFrame)
{
    m_nInSampleRate=nInSampleRate;
    m_nOutSampleRate=nOutSampleRate;
    m_nSamplesPerFrame=nSamplesPerFrame;

    m_factor = nOutSampleRate * 1.0 / nInSampleRate;
    m_nSamplesPerFrame = nSamplesPerFrame;

    if (m_high_quality)
    {
        unsigned size;

        if (m_large_filter)
            m_xoff = (LARGE_FILTER_NMULT + 1) / 2.0 * MAX(1.0,1.0/m_factor);
        else
            m_xoff = (SMALL_FILTER_NMULT + 1) / 2.0 * MAX(1.0,1.0/m_factor);

        size = (m_nSamplesPerFrame + 2*m_xoff) * sizeof(short);
        m_buffer = (short*)malloc(size);
        if (0 == m_buffer)
        {
            return -1;
        }
        memset(m_buffer,0,size);
    }
    else
    {
        m_xoff = 0;
    }
    return 0;
}

void Resample::Uninit()
{
    if (m_buffer)
    {
        free(m_buffer);
        m_buffer = 0;
    }
}

int Resample::Process(short *pInSamples,short *pOutSamples,int &nSamplesPerFrame)
{
    if(m_high_quality)
    {
        short *dst_buf;
        const short *src_buf;
        /* Okay chaps, here's how we do resampling.
         *
         * The original resample algorithm requires xoff samples *before* the
         * input buffer as history, and another xoff samples *after* the
         * end of the input buffer as lookahead. Since application can only
         * supply framesize buffer on each run, PJMEDIA needs to arrange the
         * buffer to meet these requirements.
         *
         * So here comes the trick.
         *
         * First of all, because of the history and lookahead requirement,
         * m_buffer need to accomodate framesize+2*xoff samples in its
         * buffer. This is done when the buffer is created.
         *
         * On the first run, the input frame (supplied by application) is
         * copied to m_buffer at 2*xoff position. The first 2*xoff
         * samples are initially zeroed (in the initialization). The resample
         * algorithm then invoked at m_buffer+xoff ONLY, thus giving
         * it one xoff at the beginning as zero, and one xoff at the end
         * as the end of the original input. The resample algorithm will see
         * that the first xoff samples in the input as zero.
         *
         * So here's the layout of m_buffer on the first run.
         *
         * run 0
         *     +------+------+--------------+
         *     | 0000 | 0000 |  frame0...   |
         *     +------+------+--------------+
         *     ^      ^      ^              ^
         *     0    xoff  2*xoff       size+2*xoff
         *
         * (Note again: resample algorithm is called at m_buffer+xoff)
         *
         * At the end of the run, 2*xoff samples from the end of
         * m_buffer are copied to the beginning of m_buffer.
         * The first xoff part of this will be used as history for the next
         * run, and the second xoff part of this is actually the start of
         * resampling for the next run.
         *
         * And the first run completes, the function returns.
         *
         *
         * On the next run, the input frame supplied by application is again
         * copied at 2*xoff position in the m_buffer, and the
         * resample algorithm is again invoked at m_buffer+xoff
         * position. So effectively, the resample algorithm will start its
         * operation on the last xoff from the previous frame, and gets the
         * history from the last 2*xoff of the previous frame, and the look-
         * ahead from the last xoff of current frame.
         *
         * So on this run, the buffer layout is:
         *
         * run 1
         *     +------+------+--------------+
         *     | frm0 | frm0 |  frame1...   |
         *     +------+------+--------------+
         *     ^      ^      ^              ^
         *     0    xoff  2*xoff       size+2*xoff
         *
         * As you can see from above diagram, the resampling algorithm is
         * actually called from the last xoff part of previous frame (frm0).
         *
         * And so on the process continues for the next frame, and the next,
         * and the next, ...
         *
         */
        dst_buf = m_buffer + m_xoff*2;
        memcpy(dst_buf,pInSamples,m_nSamplesPerFrame << 1);

        if (m_factor >= 1)
        {
            if (m_large_filter)
            {
                nSamplesPerFrame = SrcUp(m_buffer + m_xoff,pOutSamples,m_factor,m_nSamplesPerFrame,
                                        LARGE_FILTER_NWING,LARGE_FILTER_SCALE,
                                        LARGE_FILTER_IMP,LARGE_FILTER_IMPD,1);
            }
            else
            {
                nSamplesPerFrame = SrcUp(m_buffer + m_xoff,pOutSamples,
                                        m_factor,m_nSamplesPerFrame,SMALL_FILTER_NWING,
                                        SMALL_FILTER_SCALE,SMALL_FILTER_IMP,SMALL_FILTER_IMPD,1);
            }
        }
        else
        {
            if (m_large_filter)
            {
                nSamplesPerFrame = SrcUD(m_buffer + m_xoff,pOutSamples,m_factor, m_nSamplesPerFrame,
                                        LARGE_FILTER_NWING,LARGE_FILTER_SCALE * m_factor + 0.5,
                                        LARGE_FILTER_IMP, LARGE_FILTER_IMPD,1);
            }
            else
            {
                nSamplesPerFrame = SrcUD(m_buffer + m_xoff,pOutSamples,
                                        m_factor,m_nSamplesPerFrame,SMALL_FILTER_NWING,
                                        SMALL_FILTER_SCALE * m_factor + 0.5,SMALL_FILTER_IMP,SMALL_FILTER_IMPD,1);
            }
        }
        dst_buf = m_buffer;
        src_buf = pInSamples + m_nSamplesPerFrame - m_xoff*2;
        memcpy(dst_buf, src_buf, m_xoff * 2 * 2);
    }
    else
    {
        nSamplesPerFrame = SrcLinear(pInSamples,pOutSamples,m_factor,m_nSamplesPerFrame);
    }
    return 0;
}
