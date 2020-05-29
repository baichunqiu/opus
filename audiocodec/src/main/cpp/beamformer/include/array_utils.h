#ifndef BEAM_UTILS_H
#define BEAM_UTILS_H

#include "Platform.h"

// Cartesian Coordinate
typedef struct
{
    OsFloat x;
    OsFloat y;
    OsFloat z;
} SensorPoint;

static OsFloat Distance(SensorPoint a,SensorPoint b)
{
    return sqrt(pow(a.x-b.x,2) + pow(a.y-b.y,2) + pow(a.z-b.z,2));
}

static OsFloat GetMinimumSpacing(SensorPoint *inPos,OsInt32 inElemSize)
{
    OsFloat mic_spacing = pow(2.0,31)-1;
    OsInt32 i,j;

    for(i = 0; i < inElemSize-1; ++i)
    {
        for(j = i+1; j < inElemSize; ++j)
        {
            OsFloat distance = Distance(inPos[i],inPos[j]);
            mic_spacing = min(mic_spacing,distance);
        }
    }
    return mic_spacing;
}

static void CenteredArray(OsFloat (*ioGeometry)[3],OsInt32 inSensorNum)
{
    OsInt32 dim,i;
    for (dim = 0; dim < 3; ++dim)
    {
        OsFloat center = 0.f;
        for (i = 0; i < inSensorNum; ++i)
        {
            center += ioGeometry[i][dim];
        }
        center /= inSensorNum;
        for (i = 0; i < inSensorNum; ++i)
        {
            ioGeometry[i][dim] -= center;
        }
    }
}

static OsInt32 gcd(OsInt32 a,OsInt32 b)
{
    OsInt32 tmp;

    while (b)
    {
        tmp = a;
        a = b;
        b = tmp % b;
    }
    return a;
}

static OsFloat I0(OsFloat x)
{
    OsFloat y = x/3.75f;
    OsFloat c;

    y *= y;

    c = 1.2067492f + y * (0.2659732f + y * (0.360768e-1f + y * 0.45813e-2f));
    return 1.0f + y * (3.5156229f + y * (3.0899424f + y * c));
}

static void kaiser_c(OsInt32 n,OsFloat beta,OsFloat *w)
{
    OsInt32 half = (n + 1) / 2;
    OsFloat sum = 0.0f;
    OsInt32 i;

    for (i = 0; i <= half; ++i)
    {
        OsFloat r = (4.0f * i) / n - 1.0f;
        sum += I0(M_PI * beta * sqrt(1.0f - r * r));
        w[i] = sum;
    }
    for (i = n - 1; i >= half; --i)
    {
        w[n-i-1] = sqrtf(w[n-i-1] / sum);
        w[i] = w[n-i-1];
    }
}

#endif