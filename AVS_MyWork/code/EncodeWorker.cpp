#include "EncodeWorker.h"

EncodeWorker::EncodeWorker(int num, IRef* refptr):
    IWorker(num, refptr)
{
    m_memBufPoolPtr = m_paramPtr->m_memBufPoolPtr;

    m_outbuf.mem_malloc(m_paramPtr->m_iBlockSize/2);

}

EncodeWorker::~EncodeWorker()
{

}

void EncodeWorker::doJob()
{
    while (true)
    {
        MemBuf *ptr = m_memBufPoolPtr->getFullbuf();
        if(ptr)
        {
            if(ptr->m_len[0]==0)
            {
                delete ptr;
                break;
            }
            else
            {
                doTask(ptr);
            }
        }
        m_processcount++;
    }
}

void EncodeWorker::reorder_process(const std::vector<uint32_t> &order) 
{
    arc_stdout(level::INFO, "Reorder Process ... ");
    { // reorder seq and qual
        char * tmp = new char[m_seqElPtr->m_rawBuf.size()]; // tmp can be used by seq and qual
        std::vector<char *> tmp_start_ptr(m_seqLenElPtr->m_LenArry.size() + 1);
        tmp_start_ptr[0] = tmp;
        for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); ++i)
            tmp_start_ptr[i + 1] = tmp_start_ptr[i] + m_seqLenElPtr->m_LenArry[i];

        {
            std::memcpy(tmp, m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size());
            char *dest = m_seqElPtr->m_rawBuf.getdata();
            for (auto order_id: order) {
                std::memcpy(dest, tmp_start_ptr[order_id], m_seqLenElPtr->m_LenArry[order_id]);
                dest += m_seqLenElPtr->m_LenArry[order_id];
            }
        }

        {
            std::memcpy(tmp, m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size());
            char * dest = m_qualElPtr->m_rawBuf.getdata();
            for (auto order_id: order){
                std::memcpy(dest, tmp_start_ptr[order_id], m_seqLenElPtr->m_LenArry[order_id]);
                dest += m_seqLenElPtr->m_LenArry[order_id];
            }
        }

        delete [] tmp;
    }

    { // reorder name
        char * tmp = new char[m_nameElPtr->m_rawBuf.size()];
        std::vector<char *> tmp_start_ptr(m_namelenElPtr->m_LenArry.size() + 1);
        tmp_start_ptr[0] = tmp;
        for (size_t i = 0; i < m_namelenElPtr->m_LenArry.size(); ++i)
            tmp_start_ptr[i + 1] = tmp_start_ptr[i] + m_namelenElPtr->m_LenArry[i];
        std::memcpy(tmp, m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size());
        char * dest = m_nameElPtr->m_rawBuf.getdata();
        for (auto order_id : order){
            std::memcpy(dest, tmp_start_ptr[order_id], m_namelenElPtr->m_LenArry[order_id]);
            dest += m_namelenElPtr->m_LenArry[order_id];
        }
        delete [] tmp;
    }

    { // reorder seq_len and name_len
        if (m_seqLenElPtr->m_LenArry.size() != m_namelenElPtr->m_LenArry.size()){
            throw std::runtime_error("SeqLen number isn't equal to NameLen number");
        }
        auto * tmp = new uint32_t [m_seqLenElPtr->m_LenArry.size()];
        std::memcpy(tmp, m_seqLenElPtr->m_LenArry.getdata(), sizeof(uint32_t) * m_seqLenElPtr->m_LenArry.size());
        for (size_t i = 0; i < order.size(); ++i) {
            m_seqLenElPtr->m_LenArry[i] = tmp[order[i]];
        }
        std::memcpy(tmp, m_namelenElPtr->m_LenArry.getdata(), sizeof(uint32_t) * m_namelenElPtr->m_LenArry.size());
        for (size_t i = 0; i < order.size(); ++i){
            m_namelenElPtr->m_LenArry[i] = tmp[order[i]];
        }
        delete [] tmp;
    }
    arc_stdout(level::INFO, "Reorder Finish ... ");
}

int set_seq_compress_part(char *outptr, const char* seq_compress_part, int seq_compress_size){
    int idlen = Encap::setID(BLOCK_SEQ, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;
    std::memcpy(outptr, seq_compress_part, seq_compress_size);
    Encap::setSize(seq_compress_size, SIZENUM4, psize);
    return seq_compress_size+SIZENUM4+idlen;
}

void EncodeWorker::decomposeAlignInfo(AlignInfo &info)
{
    // int lastl = 0, remainLen = info.m_seqlen;
    // for (int i = 0; i < info.m_misnum; i++) {
    //     int2bit(info.m_cigarL[i] - lastl, getbitnum(remainLen), m_seqElPtr->m_vecCigaL);
    //     lastl = info.m_cigarL[i];
    //     remainLen = info.m_seqlen - info.m_cigarL[i];
    //     m_seqElPtr->m_vecCigaV.push_back(info.m_cigarV[i]);
    //     //printf("%d, %d  ", info.m_cigarL[i], info.m_cigarV[i]);
    // }
    // //printf("pos:%ld rev:%d mis:%d--\r\n",info.m_alignpos, info.m_bRev, info.m_misnum);

    // //m_seqElPtr->m_vecMis.push_back(info.m_misnum);
    // int2bit(info.m_misnum, getbitnum(info.m_limit), m_seqElPtr->m_vecMis);
    // m_seqElPtr->m_vecRev.push_back(info.m_bRev);
    // m_seqElPtr->m_alignCount++;
    
    //-----------------------------------------------------
    uint32_t index = 0;
    for (int i = 0; i < info.m_cigarL.size(); i++)
    {
        for (; index < info.m_cigarL[i]; index++)
        {
            m_seqElPtr->m_vecCigaL.push_back(0);
        }
        m_seqElPtr->m_vecCigaL.push_back(1);
        index++;

        m_seqElPtr->m_vecCigaV.push_back(info.m_cigarV[i]);
    }

    for (; index < info.m_seqlen; index++)
    {
        m_seqElPtr->m_vecCigaL.push_back(0);
    }

    m_seqElPtr->m_vecRev.push_back(info.m_bRev);
    m_seqElPtr->m_alignCount++;

    // int2bit(info.m_misnum, getbitnum(info.m_limit), m_seqElPtr->m_vecMis);
}

void EncodeWorker::cutHTDege(char *pseq, char *pqual, int len, AlignInfo &info)
{
    m_degeinfo.degecount = 0;
    m_degeinfo.max_dege_qual = 0;
    uint8_t tval = 0;

    int st = 0, end = len -1;
    while (hash_std_table[pseq[st]] > 3 && st<len)
    {
        m_seqElPtr->m_vecDegeCh.push_back(pseq[st]);
        m_degeinfo.degecount++;
        if(pqual[st] > m_degeinfo.max_dege_qual)
        {
            m_degeinfo.max_dege_qual = pqual[st];
        }
        st++;
    }
    info.h_degelen = st;
    while (hash_std_table[pseq[end]] > 3 && end>0)
    {
        end--;
    }
    while (st<=end)
    {
        tval = hash_std_table[pseq[st]];
        if(tval > 3)
        {
            m_seqElPtr->m_vecDegeCh.push_back(pseq[st]);
            m_degeinfo.degecount++;
            if(pqual[st] > m_degeinfo.max_dege_qual)
            {
                m_degeinfo.max_dege_qual = pqual[st];
            }
        }

        info.m_seqarry.push_back(tval);
        st++;
    }
    info.m_seqlen = info.m_seqarry.size();
    while (st<len)
    {
        m_seqElPtr->m_vecDegeCh.push_back(pseq[st]);
        m_degeinfo.degecount++;
        if(pqual[st] > m_degeinfo.max_dege_qual)
        {
            m_degeinfo.max_dege_qual = pqual[st];
        }
        st++;
    }
}

void EncodeWorker::getDegeInfo(char *pseq, int len)
{
    int degecount = 0;
    for(int i=0;i<len;i++)
    {
        if(hash_std_table[pseq[i]] > 3)
        {
            degecount++;
            m_seqElPtr->m_vecDegeCh.push_back(pseq[i]);
        }
    }
    if(degecount)
    {
        m_seqElPtr->m_vecDegeTip.push_back(1);
    }
    else
    {
        m_seqElPtr->m_vecDegeTip.push_back(0);
    }
}

void EncodeWorker::getDegeInfo(char *pseq, int len, AlignInfo &info)
{
    int degecount = 0;
    uint8_t tval = 0;
    info.m_seqlen = len;
    for(int i=0;i<len;i++)
    {
        tval = hash_std_table[pseq[i]];
        info.m_seqarry.push_back(tval);
        if(tval > 3)
        {
            degecount++;
            m_seqElPtr->m_vecDegeCh.push_back(pseq[i]);
        }
    }

    if(degecount)
    {
        m_seqElPtr->m_vecDegeTip.push_back(1);
    }
    else
    {
        m_seqElPtr->m_vecDegeTip.push_back(0);
    }
}

void EncodeWorker::getDegeInfo(char *pseq, char *pqual, int len, AlignInfo &info)
{
    m_degeinfo.degecount = 0;
    m_degeinfo.max_dege_qual = 0;
    info.m_seqlen = len;
    uint8_t tval = 0;
    for(int i=0;i<len;i++)
    {
        tval = hash_std_table[pseq[i]];
        if( tval > 3)
        {
            m_seqElPtr->m_vecDegeCh.push_back(pseq[i]);
            m_degeinfo.degecount++;
            if(pqual[i] > m_degeinfo.max_dege_qual)
            {
                m_degeinfo.max_dege_qual = pqual[i];
            }
            info.m_seqarry.push_back(tval&3);
        }
        else
        {
            info.m_seqarry.push_back(tval);
        }
    }
}

void EncodeWorker::getDegeInfo(char *pseq, char *pqual, int len)
{
    m_degeinfo.degecount = 0;
    m_degeinfo.max_dege_qual = 0;

    for(int i=0;i<len;i++)
    {
        if( hash_std_table[pseq[i]] > 3)
        {
            m_seqElPtr->m_vecDegeCh.push_back(pseq[i]);
            m_degeinfo.degecount++;
            if(pqual[i] > m_degeinfo.max_dege_qual)
            {
                m_degeinfo.max_dege_qual = pqual[i];
            }
        }
    }
}

void EncodeWorker::DegeInfoProcess(char *pseq, char *pqual, int len)
{
    if(m_degeinfo.degecount)
    {
        m_seqElPtr->m_vecDegeTip.push_back(1);
        m_seqElPtr->m_vecDegeMaxQual.push_back(m_degeinfo.max_dege_qual);

        uint32_t interval = 0;
        uint32_t cnt = 0;
        for(int i=0;i<len;i++)
        {
            if(pqual[i] <= m_degeinfo.max_dege_qual)
            {
                if(hash_std_table[pseq[i]] < 4)
                {
                    cnt++;
                    m_seqElPtr->m_vecNDegePos.push_back(interval);
                    interval = 0;
                }
                else
                {
                    interval++;
                }
            }
        }
        m_seqElPtr->m_vecNDegeCnt.push_back(cnt);
    }
    else
    {
        m_seqElPtr->m_vecDegeTip.push_back(0);
    }
}

int EncodeWorker::DoCompressData()
{
    char *outptr = m_outbuf.getdata();
    int idlen = Encap::setID(BLOCK_ENCAP, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += 8; // SIZENUM4; // block 比较大的时候,使用4字节存储已压缩的数据大小会溢出
    char *ptmp = outptr;

    int szparam, szseqlen, szname, szseq, szqual,szcombine;

    if (m_paramPtr->m_alignType) { // use reference
        szparam = DoCompressPart(m_paramProcessPtr, BLOCK_PARAM, outptr);
        outptr += szparam;
        szseqlen = DoCompressPart(m_seqlenProcessPtr, BLOCK_SEQLEN, outptr);
        outptr += szseqlen;
        szname = DoCompressPart(m_nameProcessPtr, BLOCK_NAME, outptr);
        outptr += szname;
        szseq = DoCompressPart(m_seqProcessPtr, BLOCK_SEQ, outptr);
        outptr += szseq;
        szqual = DoCompressPart(m_qualProcessPtr, BLOCK_QUAL, outptr);
        outptr += szqual;
        auto compress = m_profile->getCompress();
	    int namerawsize = m_nameElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size()*2;
	    int seqrawsize = m_seqElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size();
	    compress->recordCompress("name", namerawsize, szname);
	    compress->recordCompress("qual", seqrawsize, szqual);
	    compress->recordCompress("seq", seqrawsize, szseq+szseqlen);
	    
    } else { // ref free
        if (m_paramPtr->m_seq_compress_mode == SEQ_PRESERVE_ORDER) //保持原序输出
        {
            szparam = DoCompressPart(m_paramProcessPtr, BLOCK_PARAM, outptr);
            outptr += szparam;
            szseqlen = DoCompressPart(m_seqlenProcessPtr, BLOCK_SEQLEN, outptr);
            outptr += szseqlen;
            szname = DoCompressPart(m_nameProcessPtr, BLOCK_NAME, outptr);
            outptr += szname;

            szcombine = DoCompressPart(m_combineProcessPtr,BLOCK_COMBINE,outptr);
            outptr += szcombine;

            //szseq = DoCompressPart(m_refFreeSeqProcessPtr, BLOCK_SEQ, outptr);
            //outptr += szseq;
            //szqual = DoCompressPart(m_qualProcessPtr, BLOCK_QUAL, outptr); // 因为Qual的解压依赖Seq,所以把Qual的压缩放到Seq之后
            //outptr += szqual;
            auto compress = m_profile->getCompress();
    int namerawsize = m_nameElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size()*2;
    int seqrawsize = m_seqElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size();
    compress->recordCompress("name", namerawsize, szname);
    compress->recordCompress("combine", seqrawsize*2, szcombine+szseqlen);
    //compress->recordCompress("seq", seqrawsize, szseq+szseqlen);

        } else { //重排序
            std::vector<char> seq_compress_part(m_seqElPtr->m_rawBuf.size() / 2); // 保存碱基压缩数据的预留大小为原始数据的一半(至少要压缩一半)
            int seq_compress_size = m_refFreeSeqProcessPtr->compress(seq_compress_part.data());
            seq_compress_part.resize(seq_compress_size);
            reorder_process(m_refFreeSeqProcessPtr->order);

            szparam = DoCompressPart(m_paramProcessPtr, BLOCK_PARAM, outptr);
            outptr += szparam;
            szseqlen = DoCompressPart(m_seqlenProcessPtr, BLOCK_SEQLEN, outptr);
            outptr += szseqlen;
            szname = DoCompressPart(m_nameProcessPtr, BLOCK_NAME, outptr);
            outptr += szname;
            szseq = set_seq_compress_part(outptr, seq_compress_part.data(), seq_compress_size);
            outptr += szseq;
            szqual = DoCompressPart(m_qualProcessPtr, BLOCK_QUAL, outptr); // 因为Qual的解压依赖Seq,所以把Qual的压缩放到Seq之后
            outptr += szqual;
           auto compress = m_profile->getCompress();
    int namerawsize = m_nameElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size()*2;
    int seqrawsize = m_seqElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size();
    compress->recordCompress("name", namerawsize, szname);
    compress->recordCompress("qual", seqrawsize, szqual);
    compress->recordCompress("seq", seqrawsize, szseq+szseqlen);
        }
    }

    //auto compress = m_profile->getCompress();
    //int namerawsize = m_nameElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size()*2;
    //int seqrawsize = m_seqElPtr->m_rawBuf.size() + m_seqLenElPtr->m_LenArry.size();
    //compress->recordCompress("name", namerawsize, szname);
    //compress->recordCompress("qual", seqrawsize, szqual);
    //compress->recordCompress("seq", seqrawsize, szseq+szseqlen);
    

    int size = outptr - ptmp;
    Encap::setSize(size, 8, /* SIZENUM4,*/ psize);
    return size+8+idlen; // size+SIZENUM4+idlen;
}

int EncodeWorker::DoCompressPart(IProcess *ptr, int id, char *outptr)
{
    if (m_paramPtr->m_alignType) {
        int idlen = Encap::setID(id, outptr);
        outptr += idlen;
        char *psize = outptr;
        outptr += SIZENUM4;
        int szlen = ptr->compress(outptr);
        Encap::setSize(szlen, SIZENUM4, psize);
        return szlen + SIZENUM4 + idlen;
    } else {
        int idlen = Encap::setID(id, outptr);
        outptr += idlen;
        char *psize = outptr;
        outptr += 8;
        int szlen = ptr->compress(outptr);
        Encap::setSize(szlen, 8, psize);
        return szlen + 8 + idlen;
    }
}
