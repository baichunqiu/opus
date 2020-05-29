#ifndef COMPLEX_H
#define COMPLEX_H

#include "Platform.h"

typedef struct {
    OsFloat re;
    OsFloat im;
} complex;

static complex init_c(OsFloat re,OsFloat im)
{
    complex x;
    x.re = re;
    x.im = im;
    return x;
}

static complex add_c(complex a,complex b)
{
    complex c;
    c.re = a.re + b.re;
    c.im = a.im + b.im;
    return c;
}

static complex sub_c(complex a,complex b)
{
    complex c;
    c.re = a.re - b.re;
    c.im = a.im - b.im;
    return c;
}

static complex mul_c(complex a,complex b)
{
    complex c;
    c.re = a.re*b.re - a.im*b.im;
    c.im = a.re*b.im + a.im*b.re;
    return c;
}

static complex scale_c(complex x,OsFloat v)
{
    x.re *= v;
    x.im *= v;
    return x;
}

static complex div_c(complex n,complex d)
{
    complex c;
    OsFloat a = d.re*d.re + d.im*d.im;
    c.re = (n.re*d.re + n.im*d.im)/a;
    c.im = (n.im*d.re - n.re*d.im)/a;
    return c;
}

static complex conj_c(complex x)
{
    x.im = -x.im;
    return x;
}

static OsFloat abs_c(const complex a)
{
    return sqrt(a.re*a.re+a.im*a.im);
}

static OsFloat angle_c(const complex x)
{
    return atan2(x.im,x.re);
}

static complex sqrt_c(const complex x,OsInt32 n)
{
    OsFloat theta,A,b;
    complex c;
    
    A = abs_c(x);
    theta = angle_c(x);
    b = powl(A,1.0/n);
    c.im = b*sin(theta/n);
    c.re = b*cos(theta/n);
    return c;
}

#endif