#include "DecodePipeWorker.h"

DecodePipeWorker::DecodePipeWorker(int num,  IRef* refptr):
    DecodeWorker(num, refptr)
{
    m_memBufPoolPtr = m_paramPtr->m_memBufPoolPtr;
}

DecodePipeWorker::~DecodePipeWorker()
{
}

uint32_t DecodePipeWorker::recoverDataNoComment()
{
    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr1 = m_outbufPtr->m_bufptr[0];
    char *pbuf1 = outptr1;
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    switch (m_paramPtr->m_PipeOut)
    {
    case PIPEOUT_SE:
    case PIPEOUT_PE:
    {
        for(uint32_t i=0;i<count;i++)
        {
            len = copyReadDataNoComment(i, outptr1);
            outptr1 += len;
        }
        break;
    }
    case PIPEOUT_PE_1:
    {
        for(uint32_t i=0;i<count;i++)
        {
            if (!(i & 1)) //第一条read
            {
                len = copyReadDataNoComment(i, outptr1);
                outptr1 += len;
            }
            else
            {
                m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
                m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
                m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            }
        }
        break;
    }
    case PIPEOUT_PE_2:
    {
        for(uint32_t i=0;i<count;i++)
        {
            if (i & 1) //第二条read
            {
                len = copyReadDataNoComment(i, outptr1);
                outptr1 += len;
            } 
            else
            {
                m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
                m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
                m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            }   
        }
        break;
    }
    default:
        break;
    }

    m_outbufPtr->m_len[0] = outptr1-pbuf1;
    return m_outbufPtr->m_len[0];
}

uint32_t DecodePipeWorker::recoverDataWithComment()
{
    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr1 = m_outbufPtr->m_bufptr[0];
    char *pbuf1 = outptr1;
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    switch (m_paramPtr->m_PipeOut)
    {
    case PIPEOUT_SE:
    case PIPEOUT_PE:
    {
        for(uint32_t i=0;i<count;i++)
        {
            len = copyReadDataWithComment(i, outptr1);
            outptr1 += len;
        }
        break;
    }
    case PIPEOUT_PE_1:
    {
        for(uint32_t i=0;i<count;i++)
        {
            if (!(i & 1)) //第一条read
            {
                len = copyReadDataWithComment(i, outptr1);
                outptr1 += len;
            }
            else
            {
                m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
                m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
                m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            }
        }
        break;
    }
    case PIPEOUT_PE_2:
    {
        for(uint32_t i=0;i<count;i++)
        {
            if (i & 1) //第二条read
            {
                len = copyReadDataWithComment(i, outptr1);
                outptr1 += len;
            } 
            else
            {
                m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
                m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
                m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            }   
        }
        break;
    }
    default:
        break;
    }
    m_outbufPtr->m_len[0] = outptr1-pbuf1;
    return m_outbufPtr->m_len[0];
}

void DecodePipeWorker::outPutData(uint32_t len)
{
    m_outbufPtr->m_blockidx = m_blockParamElPtr->m_blockidx;
    m_memBufPoolPtr->addFullbuf(m_outbufPtr);
}
