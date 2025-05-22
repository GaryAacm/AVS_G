#include "EncodeSEWorker.h"


EncodeSEWorker::EncodeSEWorker(int num, IRef* refptr):
    EncodeWorker(num, refptr)
{
}

EncodeSEWorker::~EncodeSEWorker()
{
}

void EncodeSEWorker::doTask(MemBuf *ptr)
{
    auto compress = m_profile->getCompress();
    auto timer1 = m_profile->getTimer();

    memset(&m_info, 0, sizeof(BlockInfoSE));
    m_info.originalSize = ptr->m_len[0];
    m_info.originaloffset = ptr->m_offset[0];
    m_info.fileidx = ptr->m_fileidx;
    
    timer1->start("GetDataTime");
    elementClear();
    getBlockRead(ptr);
    timer1->stop();

    if (m_paramPtr->m_alignType) {
        timer1->start("AlignTime");
        doAlign();
        timer1->stop();
    }

    timer1->start("TotalCompressTime");
    m_info.compressSize = DoCompressData();
    timer1->stop();

    timer1->start("WriteTime");
    addBlockSEInfo();
    compress->recordCompress("Total", m_info.originalSize, m_info.compressSize);
    timer1->stop();
}

void EncodeSEWorker::getBlockRead(MemBuf *ptr)
{
    int idlen = 0, shlen = 0, seqlen = 0, quallen = 0;
    m_blockParamElPtr->m_onech = ptr->m_bonech;
    m_blockParamElPtr->m_blockidx = ptr->m_blockidx;
    m_sstr.clear();
    m_sstr.rdbuf()->pubsetbuf(ptr->m_bufptr[0], ptr->m_len[0]);
    char buf[1024];
    char *pseq = m_seqElPtr->m_rawBuf.getdata();
    char *pqual = m_qualElPtr->m_rawBuf.getdata();
    int seqbuflen = m_seqElPtr->m_rawBuf.capacity();
    int qualbuflen = m_qualElPtr->m_rawBuf.capacity();
    while (1)
    {
        m_sstr.getline(buf, 1024);
        if(m_sstr.eof()) break;
        idlen = m_sstr.gcount()-2;
        m_nameElPtr->m_rawBuf.append(&buf[1], idlen);
        m_namelenElPtr->m_LenArry.push_back(idlen);

        m_sstr.getline(pseq, seqbuflen);
        seqlen = m_sstr.gcount()-1;
        seqbuflen -= seqlen;
        pseq += seqlen;

        m_sstr.getline(buf, 1024);

        m_sstr.getline(pqual, qualbuflen);
        quallen = m_sstr.gcount()-1;
        qualbuflen -= quallen;
        pqual += quallen;

        if(quallen != seqlen)
        {
            arc_stdout(level::ERROR, "seqlen isnot equal quallen");
        }
        if(m_blockParamElPtr->m_longseq == 0 && seqlen > 65536)
        {
            m_blockParamElPtr->m_longseq = 1;
        }
        m_seqLenElPtr->m_LenArry.push_back(seqlen);
    }
    
    if(m_seqLenElPtr->m_LenArry.size() != m_namelenElPtr->m_LenArry.size())
    {
        arc_stdout(level::ERROR, "seqcnt isnot equal idcnt");
    }
    m_blockParamElPtr->m_readcount = m_seqLenElPtr->m_LenArry.size();
    int sz = pseq - m_seqElPtr->m_rawBuf.getdata();
    m_seqElPtr->m_rawBuf.updatesize(sz);
    m_qualElPtr->m_rawBuf.updatesize(sz);
    m_memBufPoolPtr->addEmptybuf(ptr); 
}

void EncodeSEWorker::doAlign()
{
    auto align = m_profile->getAlignment();
    int readcount = m_namelenElPtr->m_LenArry.size();
    char *pseq = m_seqElPtr->m_rawBuf.getdata();
    char *pqual = m_qualElPtr->m_rawBuf.getdata();
    uint32_t i=0;
    bool balign = false;
    if(m_paramPtr->m_alignType)
    {
        balign = true;
        int limit = readcount*0.05;
        // bool bPre = true;
        bool bPre = false;
        int refidx;
        uint64_t offset = 0;
        m_LastRefpos = 0;
        while (i<readcount)
        {
            // doAlignDivide(pseq, m_seqLenElPtr->m_LenArry[i]);
            //----------------------------------------------
            refidx = 0;
            m_alignInfo.init(m_seqLenElPtr->m_LenArry[i]);
			m_alignInfo.m_limit = m_paramPtr->m_iMaxmis*m_seqLenElPtr->m_LenArry[i]/100;
            if(m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL)
            {
                getDegeInfo(pseq, pqual, m_seqLenElPtr->m_LenArry[i], m_alignInfo);
                DegeInfoProcess(pseq, pqual, m_seqLenElPtr->m_LenArry[i]);
                pqual += m_seqLenElPtr->m_LenArry[i];
            }
            else
            {
                getDegeInfo(pseq, m_seqLenElPtr->m_LenArry[i], m_alignInfo);
            }
            m_RefPtr->doSEAlign(m_alignInfo);
            //会得到比对上的位置
            if (m_alignInfo.m_misnum != 0xffffffff) 
            {

                if (m_alignInfo.m_alignpos == m_LastRefpos)
                {
                    m_seqElPtr->m_vecPosEqual.push_back(0);
                }
                else
                {
                    m_seqElPtr->m_vecPosEqual.push_back(1);
                    m_LastRefpos = m_alignInfo.m_alignpos;
                    m_seqElPtr->m_vecBlockPos_wj.push_back(m_LastRefpos);//存储比对上的位置
                }
                decomposeAlignInfo(m_alignInfo);
                refidx = 1;
                align->updateMisOfSE(m_alignInfo.m_misnum);
            }
            m_seqElPtr->m_vecRefIdx.push_back(refidx);
            
            uint32_t len = m_seqLenElPtr->m_LenArry[i];
            if (m_alignInfo.m_misnum == 0xffffffff && len > m_paramPtr->m_MinAlignLen && len < m_paramPtr->m_MaxAlignLen)
            {
                uint32_t len1 = len / 2;
                uint32_t len2 = len - len / 2;
                char *pseq1 = pseq;
                char *pseq2 = pseq + len1;
                doAlignDivide(pseq1, len1);
                doAlignDivide(pseq2, len2);
            }
            //----------------------------------------------

            pseq += m_seqLenElPtr->m_LenArry[i];

            i++;

            if(bPre && i>limit)
            {
                bPre = false;
                if(m_seqElPtr->m_alignCount < i*0.5)
                {
                    balign = false;
                    break;
                }
            }
        }
    }

    while (i<readcount)
    {
        if(m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL)
        {
            getDegeInfo(pseq, pqual, m_seqLenElPtr->m_LenArry[i]);
            DegeInfoProcess(pseq, pqual, m_seqLenElPtr->m_LenArry[i]);
            pqual += m_seqLenElPtr->m_LenArry[i];
        }
        else
        {
            getDegeInfo(pseq, m_seqLenElPtr->m_LenArry[i]);
        }
        pseq += m_seqLenElPtr->m_LenArry[i];

        i++;
    }

    align->setBlockMapCount(m_seqElPtr->m_alignCount, readcount - m_seqElPtr->m_alignCount);
    if(balign)
    {
        align->setAlignBlockCount();
    }
}

// 再分段，递归
void EncodeSEWorker::doAlignDivide(char *pseq, uint32_t len)
{
    AlignInfo d_alignInfo;
    int refidx = 0;
    uint64_t offset = 0;
    
    d_alignInfo.init(len);
    d_alignInfo.m_limit = m_paramPtr->m_iMaxmis*len/100;
    
    getDegeInfo(pseq, len, d_alignInfo);
    m_RefPtr->doSEAlign(d_alignInfo);

    if (d_alignInfo.m_misnum != 0xffffffff) 
    {
        if (d_alignInfo.m_alignpos == m_LastRefpos)
        {
            m_seqElPtr->m_vecPosEqual.push_back(0);
        }
        else
        {
            m_seqElPtr->m_vecPosEqual.push_back(1);
            m_LastRefpos = d_alignInfo.m_alignpos;
            m_seqElPtr->m_vecBlockPos_wj.push_back(m_LastRefpos);
        }
        decomposeAlignInfo(d_alignInfo);
        refidx = 1;               
        m_seqElPtr->m_alignCount--;
    }
    m_seqElPtr->m_vecRefIdx.push_back(refidx);

    if (d_alignInfo.m_misnum == 0xffffffff && len > m_paramPtr->m_MinAlignLen)
    {
        uint32_t len1 = len / 2;
        uint32_t len2 = len - len / 2;
        char *pseq1 = pseq;
        char *pseq2 = pseq + len1;
        doAlignDivide(pseq1, len1);
        doAlignDivide(pseq2, len2);
    }
}


void EncodeSEWorker::addBlockSEInfo()
{
    std::lock_guard<std::mutex> tlock(m_mtx_BlockInfo);
    if(m_paramPtr->m_vecBlockInfoSE.empty())
    {
        m_info.compressoffset = ARCHEADLEN;
    }
    else
    {
        BlockInfoSE &last = m_paramPtr->m_vecBlockInfoSE.back();
        m_info.compressoffset = last.compressSize+last.compressoffset;
    }
    m_paramPtr->m_vecBlockInfoSE.emplace_back(m_info);
    pwrite(m_paramPtr->m_outfd[0], m_outbuf.getdata(), m_info.compressSize, m_info.compressoffset);
}






