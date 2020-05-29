
#ifndef __DATA_TYPES_H
#define __DATA_TYPES_H

typedef int int32_t;
//typedef unsigned int size_t;
//typedef int ssize_t;
typedef short int16_t;


static unsigned short
bswap_16 (unsigned short __x)
{
  return (__x >> 8) | (__x << 8);
}

static unsigned int
bswap_32 (unsigned int __x)
{
  return (bswap_16 (__x & 0xffff) << 16) | (bswap_16 (__x >> 16));
}

static unsigned long long
bswap_64 (unsigned long long __x)
{
  return (((unsigned long long) bswap_32 (__x & 0xffffffffull)) << 32) | (bswap_32 (__x >> 32));
}

#endif  // __GETOPT_H_
