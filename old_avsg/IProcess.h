#ifndef AVS_API_IPROCESS_H
#define AVS_API_IPROCESS_H

#include "Memory.h"
#include "Md5.h"
#include "Param.h"

class IProcess
{
public:
    IProcess()
    {
        m_paramPtr = Param::GetInstance();
    }
    virtual ~IProcess(){};
    virtual int compress(char *outptr) = 0;
    virtual int decompress(int inlen, char *inptr) = 0;
protected:
    Param *m_paramPtr = nullptr;
};


#endif //AVS_API_IPROCESS_H