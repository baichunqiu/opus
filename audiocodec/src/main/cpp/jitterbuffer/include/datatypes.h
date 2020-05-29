#ifndef _DATA_TYPES_H
#define _DATA_TYPES_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>


#define memmove_int16(a,b,c) (memmove((int16_t*)(a),(int16_t*)(b),sizeof(int16_t)*(c)))
#define memcpy_int16(a,b,c) (memcpy((int16_t*)(a),(int16_t*)(b),sizeof(int16_t)*(c)))
#define memset_int16(a,b,c) (memset((int16_t*)(a),(b),sizeof(int16_t)*(c)))

#ifndef PI
#define PI					3.14159265358979323846
#endif

#endif
