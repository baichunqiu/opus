//
// Created by Z-Sabine on 2018/10/22.
//

#include <cstdint>
#include <cstring>
#include <cstdio>
#include "CheckAuthenticate.h"
#include "md5.h"

#include <android/log.h>
#define TAG    "SDKHelper"
#define LOGI(...) __android_log_print(ANDROID_LOG_INFO,TAG,__VA_ARGS__)
#define LOGE(...) __android_log_print(ANDROID_LOG_ERROR,TAG,__VA_ARGS__)

CheckAuthenticate::CheckAuthenticate(uint8_t* appKey, uint32_t nKeyLen, uint8_t* appValue, uint32_t nValueLen, uint8_t* appPackage, uint32_t nPackageLen) {
    mKeyLen = nKeyLen;
    mAppKey = new uint8_t[nKeyLen+1];
    memcpy(mAppKey, appKey, nKeyLen);
    mAppKey[nKeyLen] = 0;

    mValueLen = nValueLen;
    mAppValue = new uint8_t[nValueLen+1];
    memcpy(mAppValue, appValue, nValueLen);
    mAppValue[nValueLen] = 0;

    mPackageLen = nPackageLen;
    mAppPackage = new uint8_t[nPackageLen+1];
    memcpy(mAppPackage, appPackage, nPackageLen);
    mAppPackage[nPackageLen] = 0;

    mCheckValue = CK_UNCHECKED;
}

CheckAuthenticate::~CheckAuthenticate() {
    if (mAppKey!=NULL)
        delete mAppKey;
    if (mAppValue!=NULL)
        delete mAppKey;
    if (mAppPackage!=NULL)
        delete mAppPackage;
}

bool CheckAuthenticate::checkLocak() {
    if (mCheckValue != CK_UNCHECKED)
        return static_cast<bool>(mCheckValue);

    if (mKeyLen != 32 || mValueLen != 40 || mPackageLen == 0)
        return false;

    int cs1 = checkSum(mAppKey, mKeyLen);
    char hex1[5] = {0};
    sprintf(hex1,"%x",cs1);
    int cs1_len = strlen(hex1);
    //LOGE("appKey=%s, checkSum=%d, hex=%s, size=%d", mAppKey, cs1, hex1, cs1_len);

    int cs2 = checkSum(mAppPackage, mPackageLen);
    char hex2[5] = {0};
    sprintf(hex2,"%x",cs2);
    int cs2_len = strlen(hex2);
    //LOGE("appPackage=%s, checkSum=%d, hex=%s, size=%d", mAppPackage, cs2, hex2, cs2_len);

    unsigned char av1[5] = {0}, av2[5] = {0}, md5[33] = {0};
    for (int i = 0; i < 8; ++i) {
        if (i%2 == 0)
            av1[i/2] = mAppValue[(i/2)*10];
        else
            av2[i/2] = mAppValue[(i/2)*10 + 5];
        memcpy(md5+(i*4), mAppValue+(i*5+1), 4);
    }
    //LOGE("appValue=%s, av1=%s, av2=%s, md5=%s", mAppValue, av1, av2, md5);

    if (strcmp(hex1, reinterpret_cast<const char *>(av1 + (4 - cs1_len))) == 0 && strcmp(hex2,
                                                                                         reinterpret_cast<const char *>(
                                                                                                 av2 + (4 - cs2_len))) == 0) {
        MD5_CTX context;
        MD5Init(&context);
        int md5BufferSize = cs1_len+cs2_len+mKeyLen+mPackageLen+7;
        char md5Buffer[md5BufferSize+1];
        memset(md5Buffer, 0x00, md5BufferSize+1);
        strcat(md5Buffer, "Ak");
        strcat(md5Buffer, hex1);
        strcat(md5Buffer, "=");
        strcat(md5Buffer, reinterpret_cast<const char *>(mAppKey));
        strcat(md5Buffer, ";Pa");
        strcat(md5Buffer, hex2);
        strcat(md5Buffer, "=");
        strcat(md5Buffer, reinterpret_cast<const char *>(mAppPackage));
        MD5Update(&context, reinterpret_cast<unsigned char *>(md5Buffer), md5BufferSize);
        unsigned char digest[16];
        MD5Final(digest, &context);
        unsigned char dstMd5[33] = {0};
        hex2String(reinterpret_cast<const char *>(digest), reinterpret_cast<char *>(dstMd5), 16);
        //LOGE("Md5 src=%s, dst=%s", md5Buffer, dstMd5);
        return strcmp(reinterpret_cast<const char *>(md5),
                      reinterpret_cast<const char *>(dstMd5)) == 0;
    } else
        return false;

    return true;
}

int CheckAuthenticate::checkSum(uint8_t* pData, uint32_t size) {
    int checkSums = 0;        //做累加
    if ((NULL == pData) || (0 == size))
        return checkSums;

    for (int i = 0; i < size; ++i) {
        checkSums += pData[i] & 0xff;  //b&0xff 就是将字节转为无符号类型的
    }
    return checkSums;
}

//字节流转换为十六进制字符串的另一种实现方式
void CheckAuthenticate::hex2String(const char *sSrc, char *sDest, int nSrcLen) {
    int  i;
    char szTmp[3];

    for( i = 0; i < nSrcLen; i++ )
    {
        sprintf( szTmp, "%02x", (unsigned char) sSrc[i] );
        memcpy( &sDest[i * 2], szTmp, 2 );
    }
}