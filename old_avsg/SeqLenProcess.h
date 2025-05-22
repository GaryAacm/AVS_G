#ifndef AVS_API_SEQLENPROCESS_H
#define AVS_API_SEQLENPROCESS_H

#include "IProcess.h"
#include "clr.h"
#include "simple_model.h"

class IWorker;
class SeqLenElement;
class BlockParamElement;

class SeqLenProcess:public IProcess
{
public:
    SeqLenProcess(IWorker *ptr);
    ~SeqLenProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);

private:
    int compressLen_long(char *outptr);
    int compressLen_short(char *outptr);

    int decompressLen_long(char *data);
    int decompressLen_short(char *data);
private:
    int m_last_len = 0;
    SIMPLE_MODEL<256> m_model_len1;
    SIMPLE_MODEL<256> m_model_len2;
    SIMPLE_MODEL<256> m_model_len3;
    SIMPLE_MODEL<256> m_model_len4;
    SIMPLE_MODEL<2> m_model_same_len;
    SeqLenElement *m_seqLenElPtr = nullptr;
    BlockParamElement *m_blockParamElPtr = nullptr;
    std::shared_ptr<Module> m_profile;
};



#endif