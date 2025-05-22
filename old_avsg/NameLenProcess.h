#ifndef AVS_API_NAMELENPROCESS_H
#define AVS_API_NAMELENPROCESS_H

#include "IProcess.h"

class IWorker;
class NameLenElement;

class NameLenProcess:public IProcess
{
public:
    NameLenProcess(IWorker *ptr);
    ~NameLenProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);
private:
    NameLenElement *m_nameLenElPtr = nullptr;
};


#endif //AVS_API_NAMELENPROCESS_H