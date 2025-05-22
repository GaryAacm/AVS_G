#ifndef AVS_API_EXTRACTWORKER_H
#define AVS_API_EXTRACTWORKER_H

#include "DecodeWorker.h"
#include <functional>

class ExtractWorker:public DecodeWorker
{
public:
    ExtractWorker(int num, IRef* refptr);
    ~ExtractWorker();
    
    virtual void doJob();
    virtual uint32_t recoverDataNoComment();
    virtual uint32_t recoverDataWithComment();
    virtual void outPutData(uint32_t len);
private:
    uint32_t RegexExtractSE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t RegexExtractPE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t HeadExtractSE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t HeadExtractPE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t SrandExtractSE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t SrandExtractPE(std::function<int(ExtractWorker&, int, char*)> &func);
    uint32_t GetFiles(std::function<int(ExtractWorker&, int, char*)> &func);
private:
    MemBuf *m_outbufPtr = nullptr;
    MemBufPool *m_memBufPoolPtr = nullptr;
};



#endif