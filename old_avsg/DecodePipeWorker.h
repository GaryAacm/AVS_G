/*
 * @Author: zhaozijian
 * @Date: 2021-02-08 14:21:18
 * @LastEditors: zhaozijian
 * @LastEditTime: 2021-06-04 16:43:41
 * @Description: file content
 */
#ifndef AVS_API_DECODEPIPEWORKER_H
#define AVS_API_DECODEPIPEWORKER_H

#include "DecodeWorker.h"

class DecodePipeWorker:public DecodeWorker
{
public:
    DecodePipeWorker(int num, IRef* refptr);
    ~DecodePipeWorker();

    virtual uint32_t recoverDataNoComment();
    virtual uint32_t recoverDataWithComment();
    virtual void outPutData(uint32_t len);
    virtual void outPutData_repair(uint32_t len){};
private:
    MemBuf *m_outbufPtr = nullptr;
    MemBufPool *m_memBufPoolPtr = nullptr;
};



#endif //AVS_API_DECODEPIPEWORKER_H