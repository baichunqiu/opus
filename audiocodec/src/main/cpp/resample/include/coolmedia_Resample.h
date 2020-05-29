#ifndef RESAMPLE_H
#define RESAMPLE_H

class Resample
{
public:
    Resample();
    ~Resample();
public:
    int Init(int nInSampleRate,int nOutSampleRate,int nSamplesPerFrame);
    void Uninit();
    int Process(short *pInSamples,short *pOutSamples,int &nSamplesPerFrame);
protected:
    int         m_nInSampleRate;
    int         m_nOutSampleRate;
    int         m_nSamplesPerFrame;
    double      m_factor;           /* Conversion factor = rate_out / rate_in. */
    bool        m_large_filter;     /* Large filter? */
    bool        m_high_quality;     /* Not fast? */
    unsigned    m_xoff;             /* History and lookahead size, in samples */
    short       *m_buffer;          /* Input buffer. */
};

#endif