//
// Created by Z-Sabine on 2018/10/23.
//

#ifndef SWISSSDK_MD5_H
#define SWISSSDK_MD5_H

#include "md5global.h"

/* MD5 context. */
typedef struct {
    UINT4 state[4];                                   /* state (ABCD) */
    UINT4 count[2];        /* number of bits, modulo 2^64 (lsb first) */
    unsigned char buffer[64];                         /* input buffer */
} MD5_CTX;

void MD5Init(MD5_CTX *);
void MD5Update(MD5_CTX *, unsigned char *, unsigned int);
void MD5Final(unsigned char[16], MD5_CTX *);

#endif //SWISSSDK_MD5_H
