#ifndef AVS_API_BLOCKPARAMPROCESS_H
#define AVS_API_BLOCKPARAMPROCESS_H


#include "IProcess.h"
#include "BlockParamElement.h"

enum BlockParam
{
    BLOCKPARAM_ONECH = 1,
    BLOCKPARAM_LONGSEQ,
    BLOCKPARAM_IDX,
    BLOCKPARAM_READCNT
};

class IWorker;
class BlockParamElement;

class BlockParamProcess:public IProcess
{
public:
    BlockParamProcess(IWorker *ptr);
    ~BlockParamProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);
private:
    void dodecompress(int id, int size, char *data);
private:
    BlockParamElement *m_blockParamElPtr = nullptr;
};

#endif //AVS_API_BLOCKPARAMPROCESS_H
