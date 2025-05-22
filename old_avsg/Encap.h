#ifndef AVS_API_ENCAP_H
#define AVS_API_ENCAP_H

#include "util.h"
const int SIZENUM4 = 4;
typedef struct _tagEncapInfo
{
    uint32_t len;
    uint64_t val;
}EncapInfo;


class Encap
{
public:
    Encap(/* args */);
    ~Encap();
    static void getValue(char *pbuf, EncapInfo &info);
    static int setID(uint64_t val, char *pbuf);
    static void setSize(uint64_t size, int num, char *pbuf);

    static int encapUint8(uint32_t id, uint8_t data, char *pbuf);
    static int encapUint16(uint32_t id, uint16_t data, char *pbuf);
    static int encapUint32(uint32_t id, uint32_t data, char *pbuf);
    static int encapStr(uint32_t id, string &str, char *pbuf);
    static int encapData(uint32_t id, void *pdata, int datalen, char *pbuf);
// private:
//     char *m_pbuf = nullptr;
//     //int m_len = 0;
};




#endif //AVS_API_ENCAP_H
