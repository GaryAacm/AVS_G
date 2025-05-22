#ifndef AVS_API_DECODESEWORKER_H
#define AVS_API_DECODESEWORKER_H

#include "DecodeWorker.h"

class DecodeSEWorker:public DecodeWorker
{
public:
    DecodeSEWorker(int num, IRef* refptr);
    ~DecodeSEWorker();

    virtual uint32_t recoverDataNoComment();
    virtual uint32_t recoverDataWithComment();
    virtual void outPutData(uint32_t len);
};



#endif