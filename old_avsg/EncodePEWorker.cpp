#include "EncodePEWorker.h"
#include <algorithm>


EncodePEWorker::EncodePEWorker(int num, IRef* refptr):
    EncodeWorker(num, refptr)
{
    memset(&m_insrInfo, 0, sizeof(m_insrInfo));
    memset(&m_pealign, 0, sizeof(m_pealign));
}

EncodePEWorker::~EncodePEWorker()
{
}

void EncodePEWorker::doTask(MemBuf *ptr)
{
    auto compress = m_profile->getCompress();
    auto timer1 = m_profile->getTimer();

    memset(&m_info, 0, sizeof(BlockInfoPE));
    m_info.originalSize[0] = ptr->m_len[0];
    m_info.originaloffset[0] = ptr->m_offset[0];
    m_info.originalSize[1] = ptr->m_len[1];
    m_info.originaloffset[1] = ptr->m_offset[1];

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
    addBlockPEInfo();
    compress->recordCompress("Total", m_info.originalSize[0]+m_info.originalSize[1], m_info.compressSize);
    timer1->stop();
}

void EncodePEWorker::getBlockRead(MemBuf *ptr)
{
    char *pbuf1 = ptr->m_bufptr[0];
    char *pbuf2 = ptr->m_bufptr[1];
    int idlen = 0, shlen = 0, seqlen = 0, quallen = 0;;
    int k1 = 1, k2 = 1, idx = 0;
    m_blockParamElPtr->m_onech = ptr->m_bonech;
    m_blockParamElPtr->m_blockidx = ptr->m_blockidx;
    while (idx<ptr->m_penum)
    {
        idlen = ptr->m_vecpos1[idx]- k1;
        m_nameElPtr->m_rawBuf.append(&pbuf1[k1], idlen);
        m_namelenElPtr->m_LenArry.push_back(idlen);
        k1 = ptr->m_vecpos1[idx] + 1;

        seqlen = ptr->m_vecpos1[idx+1]-k1;
        m_seqElPtr->m_rawBuf.append(&pbuf1[k1], seqlen);
        k1 = ptr->m_vecpos1[idx+2] + 1;
        if(m_blockParamElPtr->m_longseq == 0 && seqlen > 65536)
        {
            m_blockParamElPtr->m_longseq = 1;
        }

        quallen = ptr->m_vecpos1[idx+3] - k1;
        if(quallen != seqlen)
        {
            arc_stdout(level::ERROR, "1 seqlen isnot equal quallen");
        }
        m_qualElPtr->m_rawBuf.append(&pbuf1[k1], seqlen);
        m_seqLenElPtr->m_LenArry.push_back(seqlen);
        k1 = ptr->m_vecpos1[idx+3] + 2;


        idlen = ptr->m_vecpos2[idx]- k2;
        m_nameElPtr->m_rawBuf.append(&pbuf2[k2], idlen);
        m_namelenElPtr->m_LenArry.push_back(idlen);
        k2 = ptr->m_vecpos2[idx] + 1;

        seqlen = ptr->m_vecpos2[idx+1]-k2;
        m_seqElPtr->m_rawBuf.append(&pbuf2[k2], seqlen);
        k2 = ptr->m_vecpos2[idx+2] + 1;
        if(m_blockParamElPtr->m_longseq == 0 && seqlen > 65536)
        {
            m_blockParamElPtr->m_longseq = 1;
        }

        quallen = ptr->m_vecpos1[idx+3] - k2;
        if(quallen != seqlen)
        {
            arc_stdout(level::ERROR, "1 seqlen isnot equal quallen");
        }
        m_qualElPtr->m_rawBuf.append(&pbuf2[k2], quallen);
        m_seqLenElPtr->m_LenArry.push_back(seqlen);
        k2 = ptr->m_vecpos2[idx+3] + 2;
        idx += 4;
    }
    m_blockParamElPtr->m_readcount = m_seqLenElPtr->m_LenArry.size();
    m_memBufPoolPtr->addEmptybuf(ptr);
}

bool EncodePEWorker::CaclInsertSize(std::vector<int> &vec_insr)
{
    std::sort(vec_insr.begin(), vec_insr.end());
    int len = vec_insr.size();
    int rangelen = len*0.9;

    int med = 0;
    int index = 0;
    int x=0, y=0;
    int n=1;
    int dis = 0;
    int left = 0, right = 0;

 
    if(len>0)
    {
        if(len&1) //奇数
        {
            index = (len-1)/2;
            med = vec_insr[index];
            x = y = index;
        }
        else //偶数
        {
            index = (len-2)/2;
            x = index;
            y = index+1;
            med = (vec_insr[x]+vec_insr[y])/2;
        }

        while (n++)
        {
            dis = pow(2,n);
            left = med - dis;
            right = med + dis;
            while (left < vec_insr[x] && x>=0)
            {
                x--;
            }
            while (right > vec_insr[y] && y<len)
            {
                y++;
            }
            if(y-x > rangelen)
                break;
        }
        
        m_insrInfo.n = n;
        m_insrInfo.med = med;
        m_insrInfo.left = left;
        m_insrInfo.right = right;
        m_insertSize = right - left;
        m_InsrNum = n+1;
    }

    if(m_insertSize == 0)
    {
        m_insrInfo.n = 8;
        m_insrInfo.med = 256;
        m_insertSize = 512;
        m_InsrNum = 9;
    }

    m_seqElPtr->m_InsrNum = m_InsrNum;

    return true;
}

int EncodePEWorker::AlignInfoProcessPE(PEAlign &peinfo)
{
    // printf("----%ld %d %d--\r\n",peinfo.alignpos[0], peinfo.rev[0], peinfo.mis[0]);
    // printf("------%ld %d %d----\r\n",peinfo.alignpos[1], peinfo.rev[1], peinfo.mis[1]);
    int ret = -1;
    auto align = m_profile->getAlignment();
    if(peinfo.mis[0]!= 0xffffffff && peinfo.mis[1]!= 0xffffffff) //read1和read2都比对上
    {
        int index = peinfo.alignpos[0] >> m_RefPtr->m_offsetBit;
        uint64_t offset = peinfo.alignpos[0] & m_RefPtr->m_offsetSize;
        int2bit(offset, m_RefPtr->m_offsetBit, m_seqElPtr->m_vecBlockPos);

        uint64_t insr = labs(peinfo.alignpos[0] - peinfo.alignpos[1]);
        if(insr < m_insertSize) //在insert-size范围内
        {
            if(peinfo.alignpos[1] <= peinfo.alignpos[0])
            {
                m_seqElPtr->m_vecPERel.push_back(0);
            }
            else
            {
                m_seqElPtr->m_vecPERel.push_back(1);
            }
            int2bit(insr, m_InsrNum, m_seqElPtr->m_vecBlockPos);
            align->updateMisOfWithin(peinfo.mis[0]+peinfo.mis[1]);
        }
        else //在insert-size范围外
        {
            if(peinfo.alignpos[1] <= peinfo.alignpos[0])
            {
                m_seqElPtr->m_vecPERel.push_back(2);
                int2bit(peinfo.alignpos[1], getbitnum(peinfo.alignpos[0]), m_seqElPtr->m_vecBlockPos);
            }
            else
            {
                m_seqElPtr->m_vecPERel.push_back(3);
                uint64_t c1 = m_RefPtr->m_referencelen - peinfo.alignpos[0];
                int2bit(insr, getbitnum(c1), m_seqElPtr->m_vecBlockPos);
            }
            align->updateMisOutRange(peinfo.mis[0]+peinfo.mis[1]);
        }

        m_seqElPtr->m_vecRefIdx.push_back(index+1);
        m_seqElPtr->m_vecRefIdx.push_back(index+1);

    }
    else if(peinfo.mis[0]!= 0xffffffff && peinfo.mis[1]== 0xffffffff) //read1比对上，read2未比对上
    {
        int index = peinfo.alignpos[0] >> m_RefPtr->m_offsetBit;
        uint64_t offset = peinfo.alignpos[0] & m_RefPtr->m_offsetSize;
        int2bit(offset, m_RefPtr->m_offsetBit, m_seqElPtr->m_vecBlockPos);

        m_seqElPtr->m_vecRefIdx.push_back(index+1);
        m_seqElPtr->m_vecRefIdx.push_back(0);
        align->updateMisOfRead1(peinfo.mis[0]);
    }
    else if (peinfo.mis[0]== 0xffffffff && peinfo.mis[1]!= 0xffffffff) //read2比对上，read1未比对上
    {
        int index = peinfo.alignpos[1] >> m_RefPtr->m_offsetBit;
        uint64_t offset = peinfo.alignpos[1] & m_RefPtr->m_offsetSize;
        int2bit(offset, m_RefPtr->m_offsetBit, m_seqElPtr->m_vecBlockPos);

        m_seqElPtr->m_vecRefIdx.push_back(0);
        m_seqElPtr->m_vecRefIdx.push_back(index+1);
        align->updateMisOfRead2(peinfo.mis[1]);
    }
    else if(peinfo.mis[0] == 0xffffffff && peinfo.mis[1] == 0xffffffff)
    {
        m_seqElPtr->m_vecRefIdx.push_back(0);
        m_seqElPtr->m_vecRefIdx.push_back(0);
    }

    return ret;
}

void assignment(PEAlign &pealign, AlignInfo &info1, AlignInfo &info2)
{
    pealign.mis[0] = info1.m_misnum;
    pealign.rev[0] = info1.m_bRev;
    pealign.seqlen[0] = info1.m_seqlen;
    pealign.alignpos[0] = info1.m_alignpos;

    pealign.mis[1] = info2.m_misnum;
    pealign.rev[1] = info2.m_bRev;
    pealign.seqlen[1] = info2.m_seqlen;
    pealign.alignpos[1] = info2.m_alignpos;
}

void EncodePEWorker::doAlign()
{
    bool balign = false;
    auto align = m_profile->getAlignment();
    int readcount = m_seqLenElPtr->m_LenArry.size();
    uint32_t i=0;
    char *pseq1 = m_seqElPtr->m_rawBuf.getdata();
    char *pseq2 = nullptr;
    char *pqual = m_qualElPtr->m_rawBuf.getdata();
    if (m_paramPtr->m_alignType)
    {
        balign = true;
        int limit = readcount*0.05;
        bool bPre = true;
        m_vecpealign.clear();
        m_vecinsr.clear();
        m_insertSize = 0;
        while (i<readcount)
        {
            int len1 = m_seqLenElPtr->m_LenArry[i];
            int len2 = m_seqLenElPtr->m_LenArry[i+1];
            m_alignInfo1.init(len1);
            m_alignInfo2.init(len2);
            pseq2 = pseq1 + len1;
            if(m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL)
            {
                getDegeInfo(pseq1, pqual, len1, m_alignInfo1);
                DegeInfoProcess(pseq1, pqual, len1);
                pqual += len1;
                getDegeInfo(pseq2, pqual, len2, m_alignInfo2);
                DegeInfoProcess(pseq2, pqual, len2);
                pqual += len2;
            }
            else
            {
                getDegeInfo(pseq1, len1, m_alignInfo1);
                getDegeInfo(pseq2, len2, m_alignInfo2);
            }

            m_alignInfo1.m_limit = m_paramPtr->m_iMaxmis*len1/100;
            m_alignInfo2.m_limit = m_paramPtr->m_iMaxmis*len2/100;
            m_RefPtr->doPEAlign(m_alignInfo1, m_alignInfo2);
            assignment(m_pealign, m_alignInfo1, m_alignInfo2);
            if(m_alignInfo1.m_misnum != 0xffffffff)
            {
                decomposeAlignInfo(m_alignInfo1);
            }
            if(m_alignInfo2.m_misnum != 0xffffffff)
            {
                decomposeAlignInfo(m_alignInfo2);
            }

            if(m_insertSize == 0) //没有得到block对应的insRange，先把比对结果存下来
            {
                if(m_alignInfo1.m_misnum != 0xffffffff && m_alignInfo2.m_misnum != 0xffffffff)
                {
                    uint64_t insr = labs(m_alignInfo1.m_alignpos - m_alignInfo2.m_alignpos);
                    if(insr<20000) m_vecinsr.emplace_back(insr);
                }
                m_vecpealign.emplace_back(m_pealign);
            }
            else
            {
                AlignInfoProcessPE(m_pealign);
            }

            pseq1 += len1 + len2;
            i+=2;

            if(bPre && i>limit)
            {
                bPre = false;
                CaclInsertSize(m_vecinsr);
                int size = m_vecpealign.size();
                for(int j=0;j<size;j++)
                {
                    AlignInfoProcessPE(m_vecpealign[j]);
                }
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
            getDegeInfo(pseq1, pqual, m_seqLenElPtr->m_LenArry[i]);
            DegeInfoProcess(pseq1, pqual, m_seqLenElPtr->m_LenArry[i]);
            pqual += m_seqLenElPtr->m_LenArry[i];
        }
        else
        {
            getDegeInfo(pseq1, m_seqLenElPtr->m_LenArry[i]);
        }
        pseq1 += m_seqLenElPtr->m_LenArry[i];

        i++;
    }

    align->setBlockMapCount(m_seqElPtr->m_alignCount, readcount - m_seqElPtr->m_alignCount);
    if(balign)
    {
        align->setAlignBlockCount();
    }
}

void EncodePEWorker::addBlockPEInfo()
{
    std::lock_guard<std::mutex> tlock(m_mtx_BlockInfo);
    if(m_paramPtr->m_vecBlockInfoPE.empty())
    {
        m_info.compressoffset = ARCHEADLEN;
    }
    else
    {
        BlockInfoPE &last = m_paramPtr->m_vecBlockInfoPE.back();
        m_info.compressoffset = last.compressSize+last.compressoffset;
    }
    m_paramPtr->m_vecBlockInfoPE.emplace_back(m_info);
    pwrite(m_paramPtr->m_outfd[0], m_outbuf.getdata(), m_info.compressSize, m_info.compressoffset);
}