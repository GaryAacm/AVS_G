#include "ExtractWorker.h"
#include <regex.h>

ExtractWorker::ExtractWorker(int num, IRef* refptr):
    DecodeWorker(num, refptr)
{
    m_memBufPoolPtr = m_paramPtr->m_memBufPoolPtr;
}

ExtractWorker::~ExtractWorker()
{
}


uint32_t ExtractWorker::RegexExtractSE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    regex_t reg;
    regcomp(&reg, m_paramPtr->m_strPattern.c_str() , REG_EXTENDED|REG_NOSUB);
    int nmatch = 1 , ret = 0;
    regmatch_t pmatch[1]; 

    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr = m_outbufPtr->m_bufptr[0];
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        char tmpbuf[1024]={0};
        memcpy(tmpbuf, m_nameElPtr->m_bufPtr, m_namelenElPtr->m_LenArry[i]);
        
        ret = regexec(&reg, tmpbuf, nmatch, pmatch, 0);
        if(ret == 0) //匹配
        {
            len = func(*this, i, outptr);
            outptr += len;
        }
        else
        {
            m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
            m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
        }
    }
    m_outbufPtr->m_len[0] = outptr-m_outbufPtr->m_bufptr[0];
    return m_outbufPtr->m_len[0];
}

uint32_t ExtractWorker::RegexExtractPE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    regex_t reg;
    regcomp(&reg, m_paramPtr->m_strPattern.c_str() , REG_EXTENDED|REG_NOSUB);
    int nmatch = 1 , ret = 0;
    regmatch_t pmatch[1]; 

    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr1 = m_outbufPtr->m_bufptr[0];
    char *outptr2 = m_outbufPtr->m_bufptr[1];
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        char tmpbuf[1024]={0};
        memcpy(tmpbuf, m_nameElPtr->m_bufPtr, m_namelenElPtr->m_LenArry[i]);
        
        ret = regexec(&reg, tmpbuf, nmatch, pmatch, 0);
        if(ret == 0) //匹配
        {
            if (i & 1) //第二条read
            {
                len = func(*this, i, outptr2);
                outptr2 += len;
            } else //第一条read
            {
                len = func(*this, i, outptr1);
                outptr1 += len;
            }
        }
        else
        {
            m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
            m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
        }
    }
    
    int len1 = outptr1 - m_outbufPtr->m_bufptr[0];
    int len2 = outptr2 - m_outbufPtr->m_bufptr[1];
    m_outbufPtr->m_len[0] = len1;
    m_outbufPtr->m_len[1] = len2;
    return len1+len2;
}

uint32_t ExtractWorker::HeadExtractSE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr = m_outbufPtr->m_bufptr[0];
    int len = 0;
    int count = m_paramPtr->m_vecReadCnt[m_blockParamElPtr->m_blockidx];
    for(uint32_t i=0;i<count;i++)
    {
        len = func(*this, i, outptr);
        outptr += len;
    }
    m_outbufPtr->m_len[0] = outptr-m_outbufPtr->m_bufptr[0];
    return m_outbufPtr->m_len[0]; 
}

uint32_t ExtractWorker::HeadExtractPE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr1 = m_outbufPtr->m_bufptr[0];
    char *outptr2 = m_outbufPtr->m_bufptr[1];
    int len = 0;
    int count = m_paramPtr->m_vecReadCnt[m_blockParamElPtr->m_blockidx];
    for(uint32_t i=0;i<count;i++)
    {
        if (i & 1) //第二条read
        {
            len = func(*this, i, outptr2);
            outptr2 += len;
        } else //第一条read
        {
            len = func(*this, i, outptr1);
            outptr1 += len;
        }
    }
    int len1 = outptr1 - m_outbufPtr->m_bufptr[0];
    int len2 = outptr2 - m_outbufPtr->m_bufptr[1];
    m_outbufPtr->m_len[0] = len1;
    m_outbufPtr->m_len[1] = len2;
    return len1+len2;   
}

uint32_t ExtractWorker::SrandExtractSE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    int count = m_seqLenElPtr->m_LenArry.size();
    srand((unsigned)time(0)*1000);
    unordered_set<int> hash_set;
    int s_cnt = m_paramPtr->m_vecReadCnt[m_blockParamElPtr->m_blockidx];
    while (true)
    {
        hash_set.insert(rand()%count);
        if(hash_set.size() == s_cnt)
        {
            break;
        }
    }
    

    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr = m_outbufPtr->m_bufptr[0];
    int len = 0;
    int kkk = 0;
    for(uint32_t i=0;i<count;i++)
    {
        if(hash_set.find(i) != hash_set.end())
        {
            len = func(*this, i, outptr);
            outptr += len;
            kkk++;
        }
        else
        {
            m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
            m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
        }
    }
    printf("size=%d bidx=%d kkk=%d\n", hash_set.size(), m_blockParamElPtr->m_blockidx, kkk);
    m_outbufPtr->m_len[0] = outptr-m_outbufPtr->m_bufptr[0];
    return m_outbufPtr->m_len[0];
}

uint32_t ExtractWorker::SrandExtractPE(std::function<int(ExtractWorker&, int, char*)> &func)
{
    int count = m_seqLenElPtr->m_LenArry.size();
    srand(time(0));
    unordered_set<int> hash_set;
    int s_cnt = m_paramPtr->m_vecReadCnt[m_blockParamElPtr->m_blockidx];
    while (true)
    {
        hash_set.insert(rand()%count);
        if(hash_set.size() == s_cnt)
        {
            break;
        }
    }

    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr1 = m_outbufPtr->m_bufptr[0];
    char *outptr2 = m_outbufPtr->m_bufptr[1];
    int len = 0;
    for(uint32_t i=0;i<count;i++)
    {        
        if(hash_set.find(i) != hash_set.end())
        {
            if (i & 1) //第二条read
            {
                len = func(*this, i, outptr2);
                outptr2 += len;
            } else //第一条read
            {
                len = func(*this, i, outptr1);
                outptr1 += len;
            }
        }
        else
        {
            m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
            m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
            m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
        }
    }
    
    int len1 = outptr1 - m_outbufPtr->m_bufptr[0];
    int len2 = outptr2 - m_outbufPtr->m_bufptr[1];
    m_outbufPtr->m_len[0] = len1;
    m_outbufPtr->m_len[1] = len2;
    return len1+len2;
}

uint32_t ExtractWorker::GetFiles(std::function<int(ExtractWorker&, int, char*)> &func)
{
    m_outbufPtr = m_memBufPoolPtr->getEmptybuf();
    char *outptr = m_outbufPtr->m_bufptr[0];
    int len = 0;
    int count = m_seqLenElPtr->m_LenArry.size();
    for(uint32_t i=0;i<count;i++)
    {
        len = func(*this, i, outptr);
        outptr += len;
    }
    m_outbufPtr->m_len[0] = outptr-m_outbufPtr->m_bufptr[0];
    return m_outbufPtr->m_len[0];
}

uint32_t ExtractWorker::recoverDataNoComment()
{
    std::function<int(ExtractWorker&, int, char*)> func = &ExtractWorker::copyReadDataNoComment;
    switch (m_paramPtr->m_extracttype)
    {
    case EXTRACT_MATCH:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return RegexExtractPE(func);
        }
        else
        {
            return RegexExtractSE(func);
        }
        break;
    case EXTRACT_HEAD:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return HeadExtractPE(func);
        }
        else
        {
            return HeadExtractSE(func);
        }
        break;
    case EXTRACT_SRAND:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return SrandExtractPE(func);
        }
        else
        {
            return SrandExtractSE(func);
        }
        break;
    case EXTRACT_GETFILES:
        return GetFiles(func);
        break;
    default:
        break;
    }
    return 0;
}


uint32_t ExtractWorker::recoverDataWithComment()
{
    std::function<int(ExtractWorker&, int, char*)> func = &ExtractWorker::copyReadDataWithComment;
    switch (m_paramPtr->m_extracttype)
    {
    case EXTRACT_MATCH:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return RegexExtractPE(func);
        }
        else
        {
            return RegexExtractSE(func);
        }
        break;
    case EXTRACT_HEAD:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return HeadExtractPE(func);
        }
        else
        {
            return HeadExtractSE(func);
        }
        break;
    case EXTRACT_SRAND:
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            return SrandExtractPE(func);
        }
        else
        {
            return SrandExtractSE(func);
        }
        break;
    case EXTRACT_GETFILES:
        return GetFiles(func);
        break;
    default:
        break;
    }
    return 0;
}

void ExtractWorker::outPutData(uint32_t len)
{
    m_outbufPtr->m_fileidx = m_readfileidx;
    m_outbufPtr->m_blockidx = m_blockParamElPtr->m_blockidx;
    m_memBufPoolPtr->addFullbuf(m_outbufPtr);
}

void ExtractWorker::doJob()
{
    auto timer1 = m_profile->getTimer();

    int count = m_paramPtr->m_vecBlockMeta.size();
    if(m_paramPtr->m_extracttype == EXTRACT_HEAD ||
       m_paramPtr->m_extracttype == EXTRACT_SRAND)
    {
        count = m_paramPtr->m_vecReadCnt.size();
    }
    for(uint32_t i=0;i<count;i++)
    {
        if(i % m_paramPtr->m_iThreadnum != m_num)
        {
            continue;
        }

        timer1->start("ReadTime");
        elementClear();

        readData(i);
        timer1->stop();

        timer1->start("DecodeTime");
        decodeData();
        timer1->stop();

        timer1->start("RecoverTime");
        uint32_t len = 0;
        if(m_blockParamElPtr->m_onech)
        {
            len = recoverDataNoComment();
        }
        else
        {
            len = recoverDataWithComment();
        }
        timer1->stop();

        timer1->start("WriteTime");
        outPutData(len); 
        timer1->stop();
        m_processcount++;
    }

    m_file.close();
}