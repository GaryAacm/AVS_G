#ifndef AVS_API_ENCODESEWORKER_H
#define AVS_API_ENCODESEWORKER_H

#include "EncodeWorker.h"
#include <sstream>

class EncodeSEWorker:public EncodeWorker
{
public:
    EncodeSEWorker(int num, IRef* refptr);
    ~EncodeSEWorker();
    virtual void doTask(MemBuf *ptr);
    void getBlockRead(MemBuf *ptr);
    void doAlign();

    void doAlignDivide(char *pseq, uint32_t len);

    void addBlockSEInfo();

private:
    BlockInfoSE m_info;
    AlignInfo m_alignInfo;
    std::stringstream m_sstr;
};


#endif