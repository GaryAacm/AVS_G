#include "SeqLenProcess.h"
#include "IWorker.h"


SeqLenProcess::SeqLenProcess(IWorker *ptr)
{
    m_seqLenElPtr = ptr->getseqLenElPtr();
    m_blockParamElPtr = ptr->getblockparamPtr();
    m_profile = ptr->m_profile;
}

SeqLenProcess::~SeqLenProcess()
{
}

int SeqLenProcess::compress(char *outptr)
{
    auto timer1 = m_profile->getTimer();
    timer1->start("SeqLenTime");
    int sz = 0;
    m_last_len = 0;
    if(m_blockParamElPtr->m_longseq)
    {
        sz = compressLen_long(outptr);
    }
    else
    {
        sz = compressLen_short(outptr);
    }
    timer1->stop();
    return sz;
}

int SeqLenProcess::compressLen_long(char *outptr)
{
    Memory<uint32_t> &lenarry = m_seqLenElPtr->m_LenArry;
    RangeCoder rc;
    rc.output(outptr);
    rc.StartEncode();

    m_model_len1.init();
    m_model_len2.init();
    m_model_len3.init();
    m_model_len4.init();
    m_model_same_len.init();

    for (uint32_t i = 0; i < lenarry.size(); i++) {
        if (lenarry[i] != m_last_len) {
            m_model_same_len.encodeSymbol(&rc, 0);
            m_model_len1.encodeSymbol(&rc, lenarry[i] & 0xff);
            m_model_len2.encodeSymbol(&rc, (lenarry[i] >> 8) & 0xff);
            m_model_len3.encodeSymbol(&rc, (lenarry[i] >> 16) & 0xff);
            m_model_len4.encodeSymbol(&rc, (lenarry[i] >> 24) & 0xff);
            m_last_len = lenarry[i];
        } else {
            m_model_same_len.encodeSymbol(&rc, 1);
        }
    }

    rc.FinishEncode();
    return rc.size_out();
}

int SeqLenProcess::compressLen_short(char *outptr)
{
    Memory<uint32_t> &lenarry = m_seqLenElPtr->m_LenArry;
    RangeCoder rc;
    rc.output(outptr);
    rc.StartEncode();

    m_model_len1.init();
    m_model_len2.init();
    m_model_same_len.init();

    for (uint32_t i = 0; i < lenarry.size(); i++) {
        if (lenarry[i] != m_last_len) {
            m_model_same_len.encodeSymbol(&rc, 0);
            m_model_len1.encodeSymbol(&rc, lenarry[i] & 0xff);
            m_model_len2.encodeSymbol(&rc, (lenarry[i] >> 8) & 0xff);
            m_last_len = lenarry[i];
        } else {
            m_model_same_len.encodeSymbol(&rc, 1);
        }
    }

    rc.FinishEncode();
    return rc.size_out();
}

int SeqLenProcess::decompress(int inlen, char *inptr)
{
    m_last_len = 0;
    if(m_blockParamElPtr->m_longseq)
    {
        return decompressLen_long(inptr);
    }
    else
    {
        return decompressLen_short(inptr);
    }
}

int SeqLenProcess::decompressLen_long(char *data)
{
    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    m_model_len1.init();
    m_model_len2.init();
    m_model_len3.init();
    m_model_len4.init();
    m_model_same_len.init();

    int l1=0,l2=0,l3=0,l4=0;
    for (int i = 0; i < m_blockParamElPtr->m_readcount; i++) {
        if (m_model_same_len.decodeSymbol(&rc)) {
            m_seqLenElPtr->m_LenArry.push_back(m_last_len);
        } else {
            l1 = m_model_len1.decodeSymbol(&rc);
            l2 = m_model_len2.decodeSymbol(&rc);
            l3 = m_model_len3.decodeSymbol(&rc);
            l4 = m_model_len4.decodeSymbol(&rc);
            m_last_len = l1 + (l2 << 8) + (l3<<16) + (l4<<24);
            m_seqLenElPtr->m_LenArry.push_back(m_last_len);
        }
    }
    rc.FinishDecode();
    return 0;
}

int SeqLenProcess::decompressLen_short(char *data)
{
    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    m_model_len1.init();
    m_model_len2.init();
    m_model_same_len.init();

    int l1=0,l2=0,l3=0,l4=0;
    for (int i = 0; i < m_blockParamElPtr->m_readcount; i++) {
        if (m_model_same_len.decodeSymbol(&rc)) {
            m_seqLenElPtr->m_LenArry.push_back(m_last_len);
        } else {
            l1 = m_model_len1.decodeSymbol(&rc);
            l2 = m_model_len2.decodeSymbol(&rc);
            m_last_len = l1 + (l2 << 8);
            m_seqLenElPtr->m_LenArry.push_back(m_last_len);
        }
    }
    rc.FinishDecode();
    return 0;
}

