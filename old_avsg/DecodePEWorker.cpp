#include "DecodePEWorker.h"

DecodePEWorker::DecodePEWorker(int num, IRef* refptr):
    DecodeWorker(num, refptr)
{

    m_outptr1 = m_inbuf.getdata();
    m_outptr2 = m_inbuf.getdata()+m_paramPtr->m_iBlockSize/2;
}

DecodePEWorker::~DecodePEWorker()
{
}

uint32_t DecodePEWorker::recoverDataNoComment()
{
    char *outptr1 = m_outptr1;
    char *outptr2 = m_outptr2;
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        if (i & 1) //第二条read
        {
            len = copyReadDataNoComment(i, outptr2);
            outptr2 += len;
        } else //第一条read
        {
            len = copyReadDataNoComment(i, outptr1);
            outptr1 += len;
        }
    }
    int len1 = outptr1 - m_outptr1;
    int len2 = outptr2 - m_outptr2;
    return len1+len2;
}

uint32_t DecodePEWorker::recoverDataWithComment()
{    
    char *outptr1 = m_outptr1;
    char *outptr2 = m_outptr2;
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        if (i & 1) //第二条read
        {
            len = copyReadDataWithComment(i, outptr2);
            outptr2 += len;
        } else //第一条read
        {
            len = copyReadDataWithComment(i, outptr1);
            outptr1 += len;
        }
    }
    int len1 = outptr1 - m_outptr1;
    int len2 = outptr2 - m_outptr2;
    return len1+len2;
}

void DecodePEWorker::outPutData(uint32_t len)
{
    uint32_t orglen = m_blocksize[0]+m_blocksize[1];
    if(len == orglen)
    {
        pwrite(m_paramPtr->m_outfd[0], m_outptr1, m_blocksize[0], m_blockoffset[0]);
        pwrite(m_paramPtr->m_outfd[1], m_outptr2, m_blocksize[1], m_blockoffset[1]);
    }
    else
    {
        arc_stdout(level::ERROR, "decode len %ld is not equal encode len %ld", len, orglen);
    }
}