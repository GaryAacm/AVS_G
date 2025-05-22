#ifndef AVS_API_DECODEPEWORKER_H
#define AVS_API_DECODEPEWORKER_H

#include "DecodeWorker.h"

class DecodePEWorker:public DecodeWorker
{
public:
    DecodePEWorker(int num, IRef* refptr);
    ~DecodePEWorker();

    virtual uint32_t recoverDataNoComment();
    virtual uint32_t recoverDataWithComment() ;
    virtual void outPutData(uint32_t len);
private:
    char *m_outptr1 = nullptr;
    char *m_outptr2 = nullptr;
};




#endif //AVS_API_DECODEPEWORKER_H

