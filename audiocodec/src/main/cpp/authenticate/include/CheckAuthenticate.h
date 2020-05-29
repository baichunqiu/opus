//
// Created by Z-Sabine on 2018/10/22.
//

#ifndef SWISSSDK_CHECKAUTHENTICATE_H
#define SWISSSDK_CHECKAUTHENTICATE_H

#define CK_INVALID 0
#define CK_VALID 1
#define CK_UNCHECKED 2

class CheckAuthenticate {
public:

    CheckAuthenticate(uint8_t* appKey, uint32_t nKeyLen, uint8_t* appValue, uint32_t nValueLen, uint8_t* appPackage, uint32_t nPackageLen);
    ~CheckAuthenticate();

    bool checkLocak();

private:
    uint8_t* mAppKey;
    uint32_t mKeyLen;
    uint8_t* mAppValue;
    uint32_t mValueLen;
    uint8_t* mAppPackage;
    uint32_t mPackageLen;
    int mCheckValue;

private:
    int checkSum(uint8_t* pData, uint32_t len);
    void hex2String(const char *sSrc, char *sDest, int nSrcLen);
};


#endif //SWISSSDK_CHECKAUTHENTICATE_H
