#include "Encap.h"

Encap::Encap(/* args */)
{
}

Encap::~Encap()
{
}

void Encap::getValue(char *pbuf, EncapInfo &info)
{
    info.val = 0;
    uint8_t tmp = 0;
    for (int i = 7; i >= 0; i--)
    {
        tmp = 1<<i;
        if(pbuf[0] & tmp)
        {
            info.len = 8-i;
            tmp = pbuf[0] & (tmp - 1);
            info.val = tmp << (info.len-1)*8;
            break;
        }
    }

    int idx = 1;
    uint64_t tv;
    for (int i = info.len-2; i >= 0; i--,idx++)
    {
        tv = (uint8_t)pbuf[idx];
        info.val += tv<<(i*8);
    }
}

int Encap::setID(uint64_t val, char *pbuf)
{
    uint64_t tmp = 0;
    int i=1;
    for(;i<=8;i++)
    {
        tmp = (1<<(i*7))-2;
        if(val <= tmp)
        {
            val = val | (tmp+2);
            break;
        }
    }

    for (int j = 0; j < i; j++)
    {
        pbuf[j] = (val>>(i-j-1)*8) & 0xff;
    }
    
    return i;
}

/**
 * @brief  
 * @note   
 * @param  size: 后面数据的size
 * @param  num: size自身存储需要的空间
 * @param  pbuf: 存储size需要的buf
 * @retval None
 */
void Encap::setSize(uint64_t size, int num, char *pbuf)
{
    size |= (uint64_t)1<<(num*7);
    for (int i = 0; i < num; i++)
    {
        pbuf[i] = (size>>(num-i-1)*8) & 0xff;
    }
}

int Encap::encapUint8(uint32_t id, uint8_t data, char *pbuf)
{
    char *out_p = pbuf;
    int len = setID(id, out_p);
    out_p += len;
    setSize(1, 1, out_p);
    out_p += 1;

    *out_p = data;
    out_p += 1; 

    return out_p - pbuf;
}

int Encap::encapUint16(uint32_t id, uint16_t data, char *pbuf)
{
    char *out_p = pbuf;
    int len = setID(id, out_p);
    out_p += len;
    setSize(2, 1, out_p);
    out_p += 1;

    IntTo2Ch(data, out_p);
    out_p += 2; 

    return out_p - pbuf;
}

int Encap::encapUint32(uint32_t id, uint32_t data, char *pbuf)
{
    char *out_p = pbuf;
    int len = setID(id, out_p);
    out_p += len;
    setSize(4, 1, out_p);
    out_p += 1;

    IntTo4Ch(data, out_p);
    out_p += 4; 

    return out_p - pbuf;
}

int Encap::encapStr(uint32_t id, string &str, char *pbuf)
{
    char *out_p = pbuf;
    int len = setID(id, out_p);
    out_p += len;
    setSize(str.length(), 4, out_p);
    out_p += 4;

    memcpy(out_p, str.c_str(), str.length());
    out_p += str.length(); 

    return out_p - pbuf;
}

int Encap::encapData(uint32_t id, void *pdata, int datalen, char *pbuf)
{
    char *out_p = pbuf;
    int len = setID(id, out_p);
    out_p += len;
    setSize(datalen, 4, out_p);
    out_p += 4;

    memcpy(out_p, pdata, datalen);
    out_p += datalen; 

    return out_p - pbuf;
}



