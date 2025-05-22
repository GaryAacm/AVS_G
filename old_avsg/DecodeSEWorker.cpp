#include "DecodeSEWorker.h"


DecodeSEWorker::DecodeSEWorker(int num, IRef* refptr):
    DecodeWorker(num, refptr)
{
}

DecodeSEWorker::~DecodeSEWorker()
{
}

uint32_t DecodeSEWorker::recoverDataNoComment()
{
    char *outptr = m_inbuf.getdata();
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++) //一行一行复制数据
    {
        len = copyReadDataNoComment(i, outptr);
        outptr += len;
    }

    return outptr-m_inbuf.getdata();
}

uint32_t DecodeSEWorker::recoverDataWithComment()
{
    char *outptr = m_inbuf.getdata();
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        len = copyReadDataWithComment(i, outptr);
        outptr += len;
    }

    return outptr-m_inbuf.getdata();
}

void DecodeSEWorker::outPutData(uint32_t len)
{
    if(len == m_blocksize[0])
    {
        pwrite(m_paramPtr->m_outfd[0], m_inbuf.getdata(), len, m_blockoffset[0]);
    }
    else
    {
        arc_stdout(level::ERROR, "decode len %ld is not equal encode len %ld", len, m_blocksize[0]);
    }
}