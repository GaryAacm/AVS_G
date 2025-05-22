#include "BlockParamProcess.h"
#include "IWorker.h"

BlockParamProcess::BlockParamProcess(IWorker *ptr)
{
    m_blockParamElPtr = ptr->getblockparamPtr();
}

BlockParamProcess::~BlockParamProcess()
{
}


int BlockParamProcess::compress(char *outptr)
{
    char *ptmp = outptr;
    int sz = Encap::encapUint8(BLOCKPARAM_ONECH, m_blockParamElPtr->m_onech, outptr);
    outptr += sz;
    sz = Encap::encapUint8(BLOCKPARAM_LONGSEQ, m_blockParamElPtr->m_longseq, outptr);
    outptr += sz;
    sz = Encap::encapUint32(BLOCKPARAM_IDX, m_blockParamElPtr->m_blockidx, outptr);
    outptr += sz;
    sz = Encap::encapUint32(BLOCKPARAM_READCNT, m_blockParamElPtr->m_readcount, outptr);
    outptr += sz;

    return outptr - ptmp;
}

int BlockParamProcess::decompress(int inlen, char *inptr)
{
    char *data = inptr;
    EncapInfo info;
    uint32_t readlen = 0;
    while (readlen < inlen)
    {
        Encap::getValue(data, info);
        int sid = info.val;
        data += info.len;
        readlen += info.len;
        Encap::getValue(data, info);
        data += info.len;
        readlen += info.len;

        dodecompress(sid, info.val, data); 

        data += info.val;
        readlen += info.val;
    }

    return 0;
}

void BlockParamProcess::dodecompress(int id, int size, char *data)
{
    switch (id)
    {
    case BLOCKPARAM_ONECH:
        m_blockParamElPtr->m_onech = (uint8_t)data[0];
        break;
    case BLOCKPARAM_LONGSEQ:
        m_blockParamElPtr->m_longseq = (uint8_t)data[0];
        break;
    case BLOCKPARAM_IDX:
        m_blockParamElPtr->m_blockidx = DECODE_INT((unsigned char*) data);
        break;
    case BLOCKPARAM_READCNT:
        m_blockParamElPtr->m_readcount = DECODE_INT((unsigned char*) data);
        break;
    default:
        break;
    }
}