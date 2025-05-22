#include "SeqProcess.h"
#include "IWorker.h"

#define GET_SECOND(a,b,c) a > b ? (b > c ? b : ( a > c ? c : a)) : ( a > c ? a: (b > c ? c : a))
#define FILTER_MUNS 3

SeqProcess::SeqProcess(IWorker *ptr)
{
    m_refPtr = ptr->getRefPtr();
    m_seqElPtr = ptr->getseqElPtr();
    m_seqLenElPtr = ptr->getseqLenElPtr();
    m_qualElPtr = ptr->getqualElPtr();
    m_profile = ptr->m_profile;

    m_seq_size = 1 << (2 * (7+SLEVEL));
    NS_MASK = m_seq_size - 1;
    m_model_seq = new SEQ_MODEL<uint8_t>[m_seq_size];

    if(m_paramPtr->m_actionType == ACTIONTYPE_DODECODE)
    {
        for(int i=0;i<128;i++)
        {
            m_cigar_l.push_back(-1);
            m_cigar_v.push_back(-1);
        }
        m_refseq.mem_malloc(512);
    }
    

    for(int i=0;i<64;i++)
    {
        SIMPLE_MODEL<2> *ptr = new SIMPLE_MODEL<2>;
        m_vecKmodel.emplace_back(ptr);
    }
}

SeqProcess::~SeqProcess()
{
    RELEASEARRYPTR(m_model_seq);
    for(int i=0;i<64;i++)
    {
        SIMPLE_MODEL<2> *ptr = m_vecKmodel[i];
        RELEASEPTR(ptr);
    }
}

int SeqProcess::compress(char *outptr)
{
    auto timer1 = m_profile->getTimer();
    timer1->start("SeqMd5Time");
    int szmd5 = calcMd5(outptr);
    outptr += szmd5;
    timer1->stop();

    timer1->start("SeqTime");
    int szseq = 0;
    if(m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL)
    {
        szseq = docompressSeq(outptr);
    }
    else
    {
        szseq = docompressSeq_independent(outptr);
    }
    timer1->stop();
    return szmd5+szseq;
}

int SeqProcess::calcMd5(char *outptr)
{
    int idlen = Encap::setID(BLOCKSEQ_MD5, outptr);
    outptr += idlen;
    Encap::setSize(16, 1, outptr);
    outptr += 1;

    if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &m_hashval);
        memcpy(outptr, &m_hashval, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), m_szmd5);
        memcpy(outptr, m_szmd5, 16);
    }

    return 17+idlen;
}

int SeqProcess::docompressSeq(char *outptr)
{
    char *ptmp = outptr;

    kModelInit();
    int szdegetip = compressDegeTip(outptr);
    outptr += szdegetip;
    int szdege = compressDegeCh(outptr);
    outptr += szdege;
    int szdegequal = compressDegeMaxQual(outptr);
    outptr += szdegequal;
    int szdegeposc = compressNDegeCnt(outptr);
    outptr += szdegeposc;
    int szdegepos = compressNDegePos(outptr);
    outptr += szdegepos;

    if(m_paramPtr->m_alignType)
    {
        int szorder = compressRefidx(outptr);
        outptr += szorder;
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            int szinsr = compressPEInsrNum(outptr);
            outptr += szinsr;
            int szpeRel = compressPERelation(outptr);
            outptr += szpeRel;
        }
        int szpos = compressAlignInfo_Pos(outptr);
        outptr += szpos;
        int szmis = compressAlignInfo_Mis(outptr);
        outptr += szmis;
        int szrev = compressAlignInfo_Rev(outptr);
        outptr += szrev;
        int szcigal = compressAlignInfo_CigaL(outptr);//变异的部分
        outptr += szcigal;
        int szcigav = compressAlignInfo_CigaV(outptr);//变异，匹配错误的
        outptr += szcigav;
    }

    int szunmapseq = compressUnmapSeq(outptr);
    outptr += szunmapseq;

    return outptr - ptmp;
}

void SeqProcess::encode_seq(RangeCoder *rc, char *seq, int len) 
{
    int last = 0x007616c7 & NS_MASK;

    uint8_t b = 0;
    for (int i = 0; i < len; i++) { //参考10个
        b = hash_std_table[seq[i]];

        if(b<4)
        {
            m_model_seq[last].encodeSymbol(rc, b);

            last = ((last << 2) + b) & NS_MASK;
        }
    }
}

int SeqProcess::compressNDegeCnt(char *outptr)
{
    int count = m_seqElPtr->m_vecNDegeCnt.size();
    int sz = 0;
    if(count)
    {
        int idlen = Encap::setID(BLOCKSEQ_NDEGECNT, outptr);
        outptr += idlen;
        char *psize = outptr;
        outptr += SIZENUM4;

        IntTo4Ch(count, outptr);
        outptr += 4;

        RangeCoder rc;
        rc.output(outptr);
        rc.StartEncode();

        for (int i = 0; i < count; i++) {
            kModelEncode(&rc, m_seqElPtr->m_vecNDegeCnt[i]);
        }

        rc.FinishEncode();
        sz = rc.size_out();
        sz += 4;
        Encap::setSize(sz, SIZENUM4, psize);
        sz += idlen +SIZENUM4;
    }
    return sz;
}

int SeqProcess::compressNDegePos(char *outptr)
{
    int count = m_seqElPtr->m_vecNDegePos.size();
    int sz = 0;
    if(count)
    {
        int idlen = Encap::setID(BLOCKSEQ_DEGEPOS, outptr);
        outptr += idlen;
        char *psize = outptr;
        outptr += SIZENUM4;

        IntTo4Ch(count, outptr);
        outptr += 4;

        RangeCoder rc;
        rc.output(outptr);
        rc.StartEncode();

        for (int i = 0; i < count; i++) {
            kModelEncode(&rc, m_seqElPtr->m_vecNDegePos[i]);
        }
        
        rc.FinishEncode();
        sz = rc.size_out();
        sz += 4;
        Encap::setSize(sz, SIZENUM4, psize);
        sz += idlen +SIZENUM4;
    }
    return sz;
}

int SeqProcess::compressPEInsrNum(char *outptr)
{
    int len = Encap::encapUint32(BLOCKSEQ_PEINSRNUM, m_seqElPtr->m_InsrNum, outptr);
    return len;
}

int SeqProcess::compressUnmapSeq(char *outptr) 
{
    int idlen = Encap::setID(BLOCKSEQ_UNMAPSEQ, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;

    int count = m_seqLenElPtr->m_LenArry.size();    //全部行数
    RangeCoder rc;
    rc.output(outptr);
    rc.StartEncode();

    for(int i=0;i<m_seq_size;i++)
    {
        m_model_seq[i].reset();
    }
    
    char *seq_p = m_seqElPtr->m_rawBuf.getdata();
    if(m_paramPtr->m_alignType == ALIGNTYPE_REFFREE) //mis位置
    {
        for (uint32_t i = 0; i < count; i++) 
        {
            encode_seq(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
            seq_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    else
    {
        uint32_t i = 0;
        for(; i<m_seqElPtr->m_vecRefIdx.size();i++)
        {
            if(m_seqElPtr->m_vecRefIdx[i] == 0)
            {
                encode_seq(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
            }
            seq_p += m_seqLenElPtr->m_LenArry[i];
        }
        for(; i<count;i++)
        {
            encode_seq(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
            seq_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    
    rc.FinishEncode();
    int sz = rc.size_out();
    Encap::setSize(sz, SIZENUM4, psize);

    return idlen+SIZENUM4+sz;
}

int SeqProcess::compressSimpleModel(int range, int id, int count, uint8_t *parry, char *outptr)
{
    if(count)
    {
        int idlen = Encap::setID(id, outptr);
        outptr += idlen;
        char *psize = outptr;
        outptr += SIZENUM4;

        IntTo4Ch(count, outptr);
        outptr += 4;

        RangeCoder rc;
        rc.output(outptr);
        rc.StartEncode();

        switch (range)
        {
        case 2:
        {
            SIMPLE_MODEL<2> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 3:
        {
            SIMPLE_MODEL<3> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 4:
        {
            SIMPLE_MODEL<4> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 5:
        {
            SIMPLE_MODEL<5> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 8:
        {
            SIMPLE_MODEL<8> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 9:
        {
            SIMPLE_MODEL<9> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;
        case 11:
        {
            SIMPLE_MODEL<11> model;
            uint8_t b;
            for(int i=0;i<count;i++)
            {
                b = hash_std_table[parry[i]];
                model.encodeSymbol(&rc, b - 4);
            }
        }
            break;
        case QMAX:
        {
            SIMPLE_MODEL<QMAX> model;
            uint8_t q;
            for(int i=0;i<count;i++)
            {
                q = parry[i]-'!';
                model.encodeSymbol(&rc, q);
            }
        }
            break;
        case 256:
        {
            SIMPLE_MODEL<256> model;
            for(int i=0;i<count;i++)
            {
                model.encodeSymbol(&rc, parry[i]);
            }
        }
            break;    
        default:
            {
                arc_stdout(level::ERROR, "err range=%d", range);
            }
            break;
        }


        rc.FinishEncode();
        int sz = rc.size_out();
        sz += 4;
        Encap::setSize(sz, SIZENUM4, psize);

        return sz+idlen +SIZENUM4;        
    }
    return 0;
}

int SeqProcess::compressRefidx(char *outptr)
{
    int count = m_seqElPtr->m_vecRefIdx.size();
    uint8_t *parry = m_seqElPtr->m_vecRefIdx.getdata();
    return compressSimpleModel(5, BLOCKSEQ_REFIDX, count, parry, outptr);
}


int SeqProcess::compressDegeTip(char *outptr)
{
    int count = m_seqElPtr->m_vecDegeTip.size();
    uint8_t *parry = m_seqElPtr->m_vecDegeTip.getdata();
    return compressSimpleModel(2, BLOCKSEQ_DEGETIP, count, parry, outptr);
}

int SeqProcess::compressDegeCh(char *outptr) 
{
    int count = m_seqElPtr->m_vecDegeCh.size();
    uint8_t *parry = m_seqElPtr->m_vecDegeCh.getdata();
    return compressSimpleModel(11, BLOCKSEQ_DEGECH, count, parry, outptr);
}

int SeqProcess::compressDegeMaxQual(char *outptr)
{
    int count = m_seqElPtr->m_vecDegeMaxQual.size();
    uint8_t *parry = m_seqElPtr->m_vecDegeMaxQual.getdata();
    return compressSimpleModel(QMAX, BLOCKSEQ_DEGEMAXQUAL, count, parry, outptr);
}


int SeqProcess::compressPERelation(char *outptr)
{
    int count = m_seqElPtr->m_vecPERel.size();
    uint8_t *parry = m_seqElPtr->m_vecPERel.getdata();
    return compressSimpleModel(4, BLOCKSEQ_PERELATION, count, parry, outptr);
}

int SeqProcess::compressAlignInfo_Pos(char *outptr)
{
    // int count = m_seqElPtr->m_vecBlockPos.size();
    // uint8_t *parry = m_seqElPtr->m_vecBlockPos.getdata();
    // return compressSimpleModel(2, BLOCKSEQ_MAPPOS, count, parry, outptr);

    int count1 = m_seqElPtr->m_vecPosEqual.size();
    uint8_t *parry1 = m_seqElPtr->m_vecPosEqual.getdata();

    int count2 = m_seqElPtr->m_vecBlockPos_wj.size();
    uint32_t *parry2 = m_seqElPtr->m_vecBlockPos_wj.getdata();


    int idlen = Encap::setID(BLOCKSEQ_MAPPOS, outptr);
    outptr += idlen;
    char *psize = outptr; // all size
    outptr += SIZENUM4;

    IntTo4Ch(count1, outptr);
    outptr += 4;
    IntTo4Ch(count2, outptr);
    outptr += 4;
    char *count3 = outptr;
    outptr += 4;

    char *psize1 = outptr; // size1 equal part
    outptr += SIZENUM4;
    char *psize2 = outptr; // size2  Sign
    outptr += SIZENUM4;
    char *psize3 = outptr; // size3  MSB
    outptr += SIZENUM4;
    char *psize4 = outptr; // size4  LeftBit
    outptr += SIZENUM4;

    //-------------Eocoding Part 1 : PosEqual Infomation--------------
    int model_num1 = 8;
    SIMPLE_MODEL<2> * model_posequal;
    model_posequal = new SIMPLE_MODEL<2>[model_num1];
    RangeCoder rc1;
    rc1.output(outptr);
    rc1.StartEncode();

    int model_idx1 = 0;
    for (uint32_t i = 0; i < count1; i++)
    {
        model_posequal[model_idx1].encodeSymbol(&rc1, m_seqElPtr->m_vecPosEqual[i]);
        model_idx1 <<= 1;
        model_idx1 |= m_seqElPtr->m_vecPosEqual[i];
        model_idx1 &=0x7;
    }
    rc1.FinishEncode();
    int sz1 = rc1.size_out();
    outptr += sz1;

    //-------------Eocoding Part 2.1 : Compute Delta And Sign--------------
    Memory<long int> delta;
    Memory<uint8_t> Sign;
    delta.mem_malloc(count2);
    Sign.mem_malloc(count2);

    uint32_t TempPos[FILTER_MUNS];
    for (uint32_t i = 0; i < FILTER_MUNS; i++)
    {
        delta.push_back((long int)parry2[i]);
    }
    for (uint32_t i = FILTER_MUNS; i < count2; i++)
    {
        for (uint32_t j = 0; j < FILTER_MUNS; j++) TempPos[j] = parry2[i-j-1];
        sort(TempPos, TempPos + FILTER_MUNS);
        uint32_t LastPos = TempPos[(FILTER_MUNS-1)/2];

        delta.push_back((long int)parry2[i] - (long int)LastPos);
    }

    for (uint32_t i = 0; i < count2; i++)
    {
        if (delta[i] < 0)
        {
            delta[i] = -delta[i];
            Sign.push_back(1); //nagetive
        } else {
            Sign.push_back(0);  //positive
        }
    }

    //-------------Eocoding Part 2.2 : Get MSB and LeftBit-------------
    Memory<uint8_t> MSB; // range:32~0, 33 symbol total (Zero count)
    Memory<uint8_t> LeftBit; // 0 + 1 + 2 + ... + 31 + 32 = 528
    Memory<uint32_t> LB_Model_Idx;
    MSB.mem_malloc(count2);
    LeftBit.mem_malloc(count2 * 20);
    LB_Model_Idx.mem_malloc(count2 * 20);

    for (uint32_t i = 0; i < count2; i++)
    {
        int nzcount = 0;
        while (delta[i] > 0)
        {
            if (delta[i] & 1)
            {
                LeftBit.push_back(1);
            } else {
                LeftBit.push_back(0);
            }
            delta[i] >>= 1;
            nzcount++;
        }
        MSB.push_back(32 - nzcount);
    }
    for (uint32_t i = 0; i < count2; i++)
    {
        int nzcount = 32 - MSB[i];
        for (uint32_t j = nzcount*(nzcount-1)/2; j < nzcount*(nzcount+1)/2; j++)
        {
            LB_Model_Idx.push_back(j);
        }
    }

    //-------------Eocoding Part 2.3 : Encoding-------------
    SIMPLE_MODEL<2> model_Sign;
    SIMPLE_MODEL<33> * model_MSB;
    SIMPLE_MODEL<2>  * model_LeftBit;
    model_MSB = new SIMPLE_MODEL<33>[2];
    model_LeftBit = new SIMPLE_MODEL<2>[528];

    RangeCoder rc2;
    rc2.output(outptr);
    rc2.StartEncode();
    for (uint32_t i = 0; i < count2; i++)
    {
        model_Sign.encodeSymbol(&rc2, Sign[i]);
    }
    rc2.FinishEncode();
    int sz2 = rc2.size_out();
    outptr += sz2;


    RangeCoder rc3;
    rc3.output(outptr);
    rc3.StartEncode();
    for (uint32_t i = 0; i < count2; i++)
    {
        // cout << "Sign " << i << " : " << (uint32_t)Sign[i] << endl;
        model_MSB[Sign[i]].encodeSymbol(&rc3, MSB[i]);
    }
    rc3.FinishEncode();
    int sz3 = rc3.size_out();
    outptr += sz3;


    RangeCoder rc4;
    rc4.output(outptr);
    rc4.StartEncode();
    for (uint32_t i = 0; i < LeftBit.size(); i++)
    {
        model_LeftBit[LB_Model_Idx[i]].encodeSymbol(&rc4, LeftBit[i]);
    }
    rc4.FinishEncode();
    int sz4 = rc4.size_out();
    outptr += sz4;

    //-------------Eocoding Part 2.4 : Encap And Delete-------------
    IntTo4Ch(LeftBit.size(), count3);
    Encap::setSize(sz1, SIZENUM4, psize1);
    Encap::setSize(sz2, SIZENUM4, psize2);
    Encap::setSize(sz3, SIZENUM4, psize3);
    Encap::setSize(sz4, SIZENUM4, psize4);

    int sz = 4 * 3 + 4 * 4 + sz1 + sz2 + sz3 + sz4;// count * 3 + size * 4 + sz1234
    Encap::setSize(sz, SIZENUM4, psize);

    delete[] model_posequal;
    delete[] model_MSB;
    delete[] model_LeftBit;

    return sz + idlen +SIZENUM4;
}

int SeqProcess::compressAlignInfo_Mis(char *outptr)
{
    int count = m_seqElPtr->m_vecMis.size();
    uint8_t *parry = m_seqElPtr->m_vecMis.getdata();
    return compressSimpleModel(2, BLOCKSEQ_MAPMIS, count, parry, outptr);
}

int SeqProcess::compressAlignInfo_Rev(char *outptr)
{
    int count = m_seqElPtr->m_vecRev.size();
    uint8_t *parry = m_seqElPtr->m_vecRev.getdata();
    return compressSimpleModel(2, BLOCKSEQ_MAPREV, count, parry, outptr);
}

int SeqProcess::compressAlignInfo_CigaL(char *outptr)
{
    // int count = m_seqElPtr->m_vecCigaL.size();
    // uint8_t *parry = m_seqElPtr->m_vecCigaL.getdata();
    // return compressSimpleModel(2, BLOCKSEQ_MAPCIGAL, count, parry, outptr);
    
    //---------------------------------------------------------------------
    int count = m_seqElPtr->m_vecCigaL.size();

    int idlen = Encap::setID(BLOCKSEQ_MAPCIGAL, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;
    
    IntTo4Ch(count, outptr);
    outptr += 4;
    
    //--------------------
    int model_num = 1 + 11;
    SIMPLE_MODEL<2> * model_pos;
    model_pos = new SIMPLE_MODEL<2>[model_num];
    RangeCoder rc;
    rc.output(outptr);
    rc.StartEncode();

    uint32_t model_idx = 0;

    uint32_t i = 0;
    for (; i < 11 && i < count; i++)
    {
        model_idx += m_seqElPtr->m_vecCigaL[i];
        model_pos[11].encodeSymbol(&rc, m_seqElPtr->m_vecCigaL[i]);
    }
    model_idx -= m_seqElPtr->m_vecCigaL[10];

    for (; i < count; i++)
    {
        model_idx += m_seqElPtr->m_vecCigaL[i - 1];
        model_idx -= m_seqElPtr->m_vecCigaL[i - 11];
        model_pos[model_idx].encodeSymbol(&rc, m_seqElPtr->m_vecCigaL[i]);
    }
    //--------------------
    rc.FinishEncode();
    int sz = rc.size_out();
    sz += SIZENUM4;
    Encap::setSize(sz, SIZENUM4, psize);
    delete[] model_pos;

    return idlen + SIZENUM4 + sz;
}

int SeqProcess::compressAlignInfo_CigaV(char *outptr)
{
    int count = m_seqElPtr->m_vecCigaV.size();
    uint8_t *parry = m_seqElPtr->m_vecCigaV.getdata();
    return compressSimpleModel(4, BLOCKSEQ_MAPCIGAV, count, parry, outptr);
}


int SeqProcess::decompress(int inlen, char *inptr)
{
    char *data = inptr;
    EncapInfo info;
    uint32_t readlen = 0;
    kModelInit();
    while (readlen < inlen)
    {
        Encap::getValue(data, info);
        int sid = info.val;
        data += info.len;
        readlen += info.len;
        Encap::getValue(data, info);
        data += info.len;
        readlen += info.len;

        dodecompressSeq(sid, info.val, data); 

        data += info.val;
        readlen += info.val;
    }
    return 0;
}

uint64_t SeqProcess::getRealPos(bool isread2, uint64_t order)
{
    uint64_t tmppos = (order-1) << m_refPtr->m_offsetBit;
    uint64_t ret = 0;

    if(isread2)
    {
        int rel = m_seqElPtr->m_vecPERel[m_peRelidx++];
        switch (rel)
        {
        case 0: 
        {
            ret = getOffsePos(m_seqElPtr->m_InsrNum);
            m_realpos -= ret;
            break;
        }
        case 1:
        {
            ret = getOffsePos(m_seqElPtr->m_InsrNum);
            m_realpos += ret;
            break;
        }
        case 2:
        {
            m_realpos = getOffsePos(getbitnum(m_realpos));
            break;
        }
        case 3:
        {
            uint64_t c1 = m_refPtr->m_referencelen - m_realpos;
            ret = getOffsePos(getbitnum(c1));
            m_realpos += ret;
            break;
        }
        default:
            break;
        }
    }
    else
    {
        ret = getOffsePos(m_refPtr->m_offsetBit);
        m_realpos = tmppos + ret;
    }

    return m_realpos;
}

void SeqProcess::decodeAlignInfo(int seqlen, int cigar_num, uint64_t realpos)
{
    // //printf("----%ld %d %d \n", realpos, isrev, cigar_num);
    // int offset = 0;
    // int reallen = seqlen;
    // for (int i = 0; i < cigar_num; i++) {
    //     m_cigar_l[i] = getCigal(getbitnum(reallen - offset)) + offset;
    //     offset = m_cigar_l[i];
    //     m_cigar_v[i] = m_seqElPtr->m_vecCigaV[m_cigavidx++];
    //     //printf("%d--%d  ", m_cigar_l[i], m_cigar_v[i]);
    // }
    // //printf("----%ld %d %d--\n", realpos, isrev, cigar_num);

    //-------------------------------------MisPos 01 coding
    cigar_num = 0;
    int reallen = seqlen;
    for (int i = 0 ; i < reallen; i++)
    {
        if (m_cigalPtr[i] == 1)
        {
            m_cigar_l[cigar_num] = i;
            m_cigar_v[cigar_num] = m_seqElPtr->m_vecCigaV[m_cigavidx++];
            cigar_num++;
        }
    }
    m_cigalPtr += reallen;
    //-------------------------------------

    if(reallen > m_refseq.capacity())
    {
        m_refseq.mem_realloc(reallen);
    }
    m_refPtr->query(realpos, reallen, m_refseq.getdata());

    uint8_t t = 0;
    for (int i = 0; i < cigar_num; i++) {
        t = m_refseq[m_cigar_l[i]];
        m_refseq[m_cigar_l[i]] = mapBarry[t][m_cigar_v[i]];//mapvar2base(m_refseq[m_cigar_l[i]]& 0xff, m_cigar_v[i]);
    }
}

void SeqProcess::AlignInfoToSeq(bool bdege, char *pseq, char *pqual, int seqlen, int limit, uint64_t realpos)
{
    bool isrev = m_seqElPtr->m_vecRev[m_mapidx++];
    int cigar_num = getVal(m_misptr,getbitnum(limit));
    decodeAlignInfo(seqlen, cigar_num, realpos);

    if(isrev)
    {
        int j = seqlen - 1;
        for (int i = 0; i < seqlen; i++, j--) {
            pseq[i] = "ACGT"[m_refseq[j] ^ 3];
        }
    }
    else
    {
        for (int i = 0; i < seqlen; i++) {
            pseq[i] = "ACGT"[m_refseq[i]&3];
        }
    }

    if(bdege)
    {
        char max_qual = m_seqElPtr->m_vecDegeMaxQual[m_degeidx];
        uint32_t cnt = m_seqElPtr->m_vecNDegeCnt[m_degeidx];
        m_degeidx++;

        if(cnt)
        {
            uint32_t ndegepos = m_seqElPtr->m_vecNDegePos[m_ndegeposidx++]; 
            cnt--;
            uint32_t interval = 0;
            for(int i=0;i<seqlen;i++)
            {
                if(pqual[i] <= max_qual)
                {
                    if(interval == ndegepos)
                    {
                        interval = 0;
                        if(cnt)
                        {
                            ndegepos = m_seqElPtr->m_vecNDegePos[m_ndegeposidx++];
                            cnt--;
                        }
                        else
                        {
                            ndegepos = 0xffffffff;
                        }
                    }
                    else
                    {
                        interval++;
                        pseq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                    }
                }
            }
        }
        else
        {
            for(int i=0;i<seqlen;i++)
            {
                if(pqual[i] <= max_qual)
                {
                    pseq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                }
            }
        }
    }
}

uint64_t SeqProcess::getOffsePos(int count)
{
    uint64_t ret = 0;
    uint64_t tmp;
    for(int i=0;i<count;i++)
    {
        tmp = m_posPtr[i];
        ret += tmp<<i;
    }
    m_posPtr += count;
    return ret;
}

int SeqProcess::getCigal(int count)
{
    int ret = 0;
    int tmp;
    for(int i=0;i<count;i++)
    {
        tmp = m_cigalPtr[i];
        ret += tmp<<i;
    }
    m_cigalPtr += count;
    return ret;
}


int SeqProcess::dodecompressSeq(int id, uint32_t size, char* data)
{
    switch (id)
    {
    case BLOCKSEQ_MD5:
        if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
        {
            memcpy(&m_hashval, data, 16);
        }
        else if(m_paramPtr->m_checktype == DATACHECK_MD5)
        {
            memcpy(m_szmd5, data, 16);
        }
        break;
    case BLOCKSEQ_REFIDX:
        decompressRefidx(data);
        break;
    case BLOCKSEQ_DEGETIP:
        decompressDegeTip(data);
        break;
    case BLOCKSEQ_DEGECH:
        decompressDegeCh(data);
        break;
    case BLOCKSEQ_DEGEMAXQUAL:
        decompressDegeMaxQual(data);
        break;
    case BLOCKSEQ_NDEGECNT:
        decompressNDegeCnt(data);
        break;
    case BLOCKSEQ_DEGEPOS:
        decompressNDegePos(data);
        break;
    case BLOCKSEQ_PEINSRNUM:
        decompressPEInsrNum(data);
        break;
    case BLOCKSEQ_PERELATION:
        decompressPERelation(data);
        break;
    case BLOCKSEQ_MAPPOS:
        decompressAlignInfo_Pos(data);
        break;
    case BLOCKSEQ_MAPMIS:
        decompressAlignInfo_Mis(data);
        break;
    case BLOCKSEQ_MAPREV:
        decompressAlignInfo_Rev(data);
        break;
    case BLOCKSEQ_MAPCIGAL:
        decompressAlignInfo_CigaL(data);
        break;
    case BLOCKSEQ_MAPCIGAV:
        decompressAlignInfo_CigaV(data);
        break;
    case BLOCKSEQ_UNMAPSEQ:
    {
        if(m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL)
        {
            decompressUnmapSeq(data);
        }
        else
        {
            decompressUnmapSeq_independent(data);
        }
        checkData();
        break;
    }
    default:
        break;
    }

    return 0;
}

int SeqProcess::decompressSimpleModel(int range, char* data, Memory<uint8_t> &arry)
{
    int count = DECODE_INT((unsigned char *)data);
    data += 4;

    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    uint8_t ret = 0;
    switch (range)
    {
    case 2:
    {
        SIMPLE_MODEL<2> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        }
    }
        break;
    case 3:
    {
        SIMPLE_MODEL<3> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        }
    }
        break;
    case 4:
    {
        SIMPLE_MODEL<4> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        } 
    }
        break;
    case 5:
    {
        SIMPLE_MODEL<5> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        }
    }
        break;
    case 8:
    {
        SIMPLE_MODEL<8> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        }
    }    
        break;
    case 9:
    {
        SIMPLE_MODEL<9> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc);
            arry.push_back(ret);
        }
    }
        break;
    case 11:
    {
        SIMPLE_MODEL<11> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc)+4;
            arry.push_back(fqzdec[ret]);
        }
    }
        break;
    case QMAX:
    {
        SIMPLE_MODEL<QMAX> model;
        for (int i = 0; i < count; i++) 
        {
            ret = model.decodeSymbol(&rc)+'!';
            arry.push_back(ret);
        }
    }
        break;
    default:
    {
        arc_stdout(level::ERROR, "err range=%d", range);
    }
        break;
    }
    rc.FinishDecode();
    return 0;
}

int SeqProcess::decompressRefidx(char* data)
{
    return decompressSimpleModel(5, data, m_seqElPtr->m_vecRefIdx);
}

int SeqProcess::decompressDegeTip(char* data)
{
    return decompressSimpleModel(2, data, m_seqElPtr->m_vecDegeTip);
}

int SeqProcess::decompressDegeCh(char* data)
{
    return decompressSimpleModel(11, data, m_seqElPtr->m_vecDegeCh);
}

int SeqProcess::decompressDegeMaxQual(char* data)
{
    return decompressSimpleModel(QMAX, data, m_seqElPtr->m_vecDegeMaxQual);
}

int SeqProcess::decompressPEInsrNum(char* data)
{
    m_seqElPtr->m_InsrNum = DECODE_INT((unsigned char *)data);
    return 0;
}

int SeqProcess::decompressPERelation(char* data)
{
    return decompressSimpleModel(4, data, m_seqElPtr->m_vecPERel);
}

int SeqProcess::decompressAlignInfo_Pos(char* data)
{
    // return decompressSimpleModel(2, data, m_seqElPtr->m_vecBlockPos);
    //--------------------------------------------

    int count1 = DECODE_INT((unsigned char *)data);
    data += 4;
    int count2 = DECODE_INT((unsigned char *)data);
    data += 4;
    int count3 = DECODE_INT((unsigned char *)data);
    data += 4;

    EncapInfo info;
    Encap::getValue(data, info);
    data += info.len;
    int psize1 = info.val;
    Encap::getValue(data, info);
    data += info.len;
    int psize2 = info.val;
    Encap::getValue(data, info);
    data += info.len;
    int psize3 = info.val;
    Encap::getValue(data, info);
    data += info.len;
    int psize4 = info.val;

    //-------------Decoding Part 1 : PosEqual Infomation--------------
    int model_num1 = 8;
    SIMPLE_MODEL<2> * model_posequal;
    model_posequal = new SIMPLE_MODEL<2>[model_num1];

    RangeCoder rc1;
    rc1.input(data);
    rc1.StartDecode();
    uint8_t ret1 = 0;
    int model_idx1 = 0;
    for (int i = 0; i < count1; i++)
    {
        ret1 = model_posequal[model_idx1].decodeSymbol(&rc1);
        m_seqElPtr->m_vecPosEqual.push_back(ret1);
        model_idx1 = ((model_idx1 << 1) | m_seqElPtr->m_vecPosEqual[i]) & 0x7;
    }
    rc1.FinishDecode();
    data += psize1;

    //-------------Docoding Part 2.1 : Decode Sign And MSB-------------
    Memory<uint8_t> Sign;
    Memory<uint8_t> MSB; // range:32~0, 33 symbol total (Zero count)
    Sign.mem_malloc(count2);
    MSB.mem_malloc(count2);

    SIMPLE_MODEL<2> model_Sign;
    SIMPLE_MODEL<33> * model_MSB;
    model_MSB = new SIMPLE_MODEL<33>[2];

    RangeCoder rc2;
    rc2.input(data);
    rc2.StartDecode();
    uint8_t ret2 = 0;
    for (uint32_t i = 0; i < count2; i++)
    {
        ret2 = model_Sign.decodeSymbol(&rc2);
        Sign.push_back(ret2);
    }
    rc2.FinishDecode();
    data += psize2;


    RangeCoder rc3;
    rc3.input(data);
    rc3.StartDecode();
    uint8_t ret3 = 0;
    for (uint32_t i = 0; i < count2; i++)
    {
        ret3 = model_MSB[Sign[i]].decodeSymbol(&rc3);
        MSB.push_back(ret3);
    }
    rc3.FinishDecode();
    data += psize3;


    //-------------Docoding Part 2.2 : Generate LB_Model_Idx by MSB-------------
    Memory<uint32_t> LB_Model_Idx;
    LB_Model_Idx.mem_malloc(count2 * 20);

    for (uint32_t i = 0; i < count2; i++)
    {
        int nzcount = 32 - MSB[i];
        for (uint32_t j = nzcount*(nzcount-1)/2; j < nzcount*(nzcount+1)/2; j++)
        {
            LB_Model_Idx.push_back(j);
        }
    }


    //-------------Docoding Part 2.3 : Decode LeftBit-------------
    Memory<uint8_t> LeftBit; // 0 + 1 + 2 + ... + 31 + 32 = 528
    LeftBit.mem_malloc(count2 * 20);
    SIMPLE_MODEL<2>  * model_LeftBit;
    model_LeftBit = new SIMPLE_MODEL<2>[528];

    RangeCoder rc4;
    rc4.input(data);
    rc4.StartDecode();
    uint8_t ret4 = 0;
    for (uint32_t i = 0; i < count3; i++)
    {
        ret4 = model_LeftBit[LB_Model_Idx[i]].decodeSymbol(&rc4);
        LeftBit.push_back(ret4);
    }
    rc4.FinishDecode();


    //-------------Docoding Part 2.4 : Compose to Get delta and RefPos-------------
    uint32_t LB_index = 0;
    Memory<long int> delta;
    delta.mem_malloc(count2);
    for (uint32_t i = 0; i < count2; i++)
    {
        long int deltaTemp = 0;
        int nzcount = 32 - MSB[i];
        for (int j = nzcount - 1; j >=0; j--)
        {
            deltaTemp <<=1;
            deltaTemp |= LeftBit[LB_index + j];
        }
        LB_index += nzcount;
        delta.push_back(deltaTemp);
    }

    for (uint32_t i = 0; i < count2; i++)
    {
        if (Sign[i] == 1)
        {
            delta[i] = -delta[i];
        }
    }

    uint32_t TempPos[FILTER_MUNS];
    uint32_t* parry2 = m_seqElPtr->m_vecBlockPos_wj.getdata();
    for (uint32_t i = 0; i < FILTER_MUNS; i++)
    {
        m_seqElPtr->m_vecBlockPos_wj.push_back((uint32_t)delta[i]);
    }
    for (uint32_t i = FILTER_MUNS; i < count2; i++)
    {   
        for (uint32_t j = 0; j < FILTER_MUNS; j++) TempPos[j] = parry2[i-j-1];
        sort(TempPos, TempPos + FILTER_MUNS);
        long int PosFront = TempPos[(FILTER_MUNS-1)/2];

        PosFront += delta[i];
        m_seqElPtr->m_vecBlockPos_wj.push_back((uint32_t)PosFront);
        parry2 = m_seqElPtr->m_vecBlockPos_wj.getdata(); // 更新parry2，防止realloc导致内存位置改变
    }

    
    //-------------Docoding Part 3 : Trans to Bit-------------
    uint32_t PosIndex = 0; // for each uniqe refpos
    uint32_t PEIndex = 0; // total align read nums
    for (uint32_t i = 0; i < m_seqElPtr->m_vecRefIdx.size(); i++)
    {
        if (m_seqElPtr->m_vecRefIdx[i] != 0)
        {
            if (m_seqElPtr->m_vecPosEqual[PEIndex++] == 1 && PEIndex != 1)
            {
                PosIndex++;
            }
            uint32_t PosTemp = m_seqElPtr->m_vecBlockPos_wj[PosIndex];

            m_seqElPtr->m_vecRefIdx[i] += (PosTemp >> 30);
            int2bit(PosTemp & ((1 << 30) - 1), 30, m_seqElPtr->m_vecBlockPos);//针对64bit reference 待改进
        }
    }

    
    delete[] model_posequal;
    delete[] model_MSB;
    delete[] model_LeftBit;


    return 0;

    //bit2int
    // for (int i = 0; i < count2; i++)
    // {
    //     uint32_t ret3 = 0;
    //     for (int index = 32 * i; index < 32 * (i + 1); index++)
    //     {
    //         ret3 <<= 1;
    //         ret3 |= BitPlant[index];
    //     }
    //     m_seqElPtr->m_vecBlockPos_wj.push_back(ret3);
    // }


}

int SeqProcess::decompressAlignInfo_Mis(char* data)
{
    return decompressSimpleModel(2, data, m_seqElPtr->m_vecMis);
}

int SeqProcess::decompressAlignInfo_Rev(char* data)
{
    return decompressSimpleModel(2, data, m_seqElPtr->m_vecRev);
}

int SeqProcess::decompressAlignInfo_CigaL(char* data)
{
    // return decompressSimpleModel(2, data, m_seqElPtr->m_vecCigaL);
    

    int count = DECODE_INT((unsigned char *)data);
    data += 4;

    //----------------------Range Coding--------------------
    int model_num = 1 + 11;
    SIMPLE_MODEL<2> * model_pos;
    model_pos = new SIMPLE_MODEL<2>[model_num];
    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    uint8_t ret = 0;
    uint32_t model_idx = 0;
    
    int i = 0;
    for (; i < 11 && i < count; i++)
    {
        ret = model_pos[11].decodeSymbol(&rc);
        m_seqElPtr->m_vecCigaL.push_back(ret);

        model_idx += m_seqElPtr->m_vecCigaL[i];
    }
    model_idx -= m_seqElPtr->m_vecCigaL[10];

    for (; i < count; i++)
    {
        model_idx += m_seqElPtr->m_vecCigaL[i - 1];
        model_idx -= m_seqElPtr->m_vecCigaL[i - 11];

        ret = model_pos[model_idx].decodeSymbol(&rc);
        m_seqElPtr->m_vecCigaL.push_back(ret);
    }
    rc.FinishDecode();
    return 0;
}

int SeqProcess::decompressAlignInfo_CigaV(char* data)
{
    return decompressSimpleModel(4, data, m_seqElPtr->m_vecCigaV);
}

int SeqProcess::decompressNDegeCnt(char *data)
{
    int count = DECODE_INT((unsigned char *)data);
    data += 4;

    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();
    uint32_t ret;
    for (int i = 0; i < count; i++) {
        ret = kModelDecode(&rc);
        m_seqElPtr->m_vecNDegeCnt.push_back(ret);
    }
    rc.FinishDecode();
    return 0;
}

int SeqProcess::decompressNDegePos(char *data)
{
    int count = DECODE_INT((unsigned char *)data);
    data += 4;

    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();
    uint32_t ret;
    for (int i = 0; i < count; i++) {
        ret = kModelDecode(&rc);
        m_seqElPtr->m_vecNDegePos.push_back(ret);
    }
    rc.FinishDecode();
    return 0;
}

int SeqProcess::decompressUnmapSeq(char* data)
{
    int count = m_seqLenElPtr->m_LenArry.size();
    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    for(int i=0;i<m_seq_size;i++){
        m_model_seq[i].reset();
    }
    
    char *seq_p = m_seqElPtr->m_rawBuf.getdata();
    char *qual_p = m_qualElPtr->m_rawBuf.getdata();
    m_degeidx = 0;
    m_degechidx = 0;
    m_ndegeposidx = 0;
    char *ptmp = seq_p;

    if(m_paramPtr->m_alignType == ALIGNTYPE_REFFREE)
    {
        for (uint32_t i = 0; i < count; i++) {
            decode_seq(m_seqElPtr->m_vecDegeTip[i], &rc, seq_p, qual_p, m_seqLenElPtr->m_LenArry[i]);
            seq_p += m_seqLenElPtr->m_LenArry[i];
            qual_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    else
    {
        uint32_t i = 0, limit;
        uint64_t realpos = 0;
        bool ret = false;
        m_posPtr = m_seqElPtr->m_vecBlockPos.getdata();
        m_cigalPtr = m_seqElPtr->m_vecCigaL.getdata();
        m_misptr = m_seqElPtr->m_vecMis.getdata();
        m_mapidx = 0;
        m_cigavidx = 0;
        m_peRelidx = 0;

        for (;i < m_seqElPtr->m_vecRefIdx.size(); i++) {
            if(m_seqElPtr->m_vecRefIdx[i])
            {
                ret = false;
                limit = m_seqLenElPtr->m_LenArry[i]*m_paramPtr->m_iMaxmis/100;
                if(m_paramPtr->m_arcType == ARCTYPE_PE && (i & 1) && m_seqElPtr->m_vecRefIdx[i-1]>0)
                {
                    ret = true;
                }
                realpos = getRealPos(ret, m_seqElPtr->m_vecRefIdx[i]);
                AlignInfoToSeq(m_seqElPtr->m_vecDegeTip[i], seq_p, qual_p, m_seqLenElPtr->m_LenArry[i], limit, realpos);
            }
            else
            {
                decode_seq(m_seqElPtr->m_vecDegeTip[i], &rc, seq_p, qual_p, m_seqLenElPtr->m_LenArry[i]);
            }
            
            seq_p += m_seqLenElPtr->m_LenArry[i];
            qual_p += m_seqLenElPtr->m_LenArry[i];
        }

        for (;i < count; i++) {
            decode_seq(m_seqElPtr->m_vecDegeTip[i], &rc, seq_p, qual_p, m_seqLenElPtr->m_LenArry[i]);
            seq_p += m_seqLenElPtr->m_LenArry[i];
            qual_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    
    rc.FinishDecode();
    int size = seq_p - ptmp;
    m_seqElPtr->m_rawBuf.updatesize(size);
    return size;
}

void SeqProcess::decode_seq(bool bdege, RangeCoder *rc, char *seq, char *pqual, int len)
{
    int last = 0x007616c7 & NS_MASK;
    uint8_t b = 0;
    if(bdege)
    {
        char max_qual = m_seqElPtr->m_vecDegeMaxQual[m_degeidx];
        uint32_t cnt = m_seqElPtr->m_vecNDegeCnt[m_degeidx];
        m_degeidx++;

        if(cnt)
        {
            uint32_t ndegepos = m_seqElPtr->m_vecNDegePos[m_ndegeposidx++]; 
            cnt--;
            uint32_t interval = 0;
            for(int i=0;i<len;i++)
            {
                if(pqual[i] <= max_qual)
                {
                    if(interval == ndegepos)
                    {
                        interval = 0;
                        if(cnt)
                        {
                            ndegepos = m_seqElPtr->m_vecNDegePos[m_ndegeposidx++];
                            cnt--;
                        }
                        else
                        {
                            ndegepos = 0xffffffff;
                        }
                        
                        b = m_model_seq[last].decodeSymbol(rc);
                        last = (last * 4 + b) & NS_MASK;

                        seq[i] = fqzdec[b];
                    }
                    else
                    {
                        interval++;
                        seq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                    }
                }
                else
                {
                    b = m_model_seq[last].decodeSymbol(rc);
                    last = (last * 4 + b) & NS_MASK;

                    seq[i] = fqzdec[b];
                }
            }
        }
        else
        {
            for(int i=0;i<len;i++)
            {
                if(pqual[i] <= max_qual)
                {
                    seq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                }
                else
                {
                    b = m_model_seq[last].decodeSymbol(rc);
                    last = (last * 4 + b) & NS_MASK;

                    seq[i] = fqzdec[b];
                }
            }
        }
    }
    else
    {
        for (int i = 0; i < len; i++) {
            b = m_model_seq[last].decodeSymbol(rc);
            last = (last * 4 + b) & NS_MASK;

            seq[i] = fqzdec[b];
        }
    }
}

void SeqProcess::kModelInit()
{
    m_lenModel.init();
    for(int i=0;i<64;i++)
    {
        m_vecKmodel[i]->init();
    }
}

void SeqProcess::kModelEncode(RangeCoder *rc, uint64_t val)
{
    int count = 0;
    while (val)
    {
        m_bitarry[count] = val & 1;
        val >>= 1;
        count++;
    }

    m_lenModel.encodeSymbol(rc, count);
    for(int i=0;i<count;i++)
    {
        m_vecKmodel[i]->encodeSymbol(rc, m_bitarry[i]);
    }
}

uint64_t SeqProcess::kModelDecode(RangeCoder *rc)
{
    int count = m_lenModel.decodeSymbol(rc);
    uint64_t ret = 0;
    uint64_t tmp;
    for(int i=0;i<count;i++)
    {
        tmp = m_vecKmodel[i]->decodeSymbol(rc);
        ret += tmp<<i;
    }

    return ret;
}

bool SeqProcess::checkData()
{
    if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        XXH128_hash_t res;
        bool ret = calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &res);
        if(!ret || XXH128_isEqual(res, m_hashval)==0)
        {
            arc_stdout(level::ERROR, "name md5 check fail");
        }
    }
    else if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        unsigned char szmd5[16] = {0};
        MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), szmd5);

        if(memcmp(m_szmd5, szmd5, 16) != 0)
        {
            // FILE *f = fopen("err.txt", "w");
            // char *ptr = m_seqElPtr->m_rawBuf.getdata();
            // for(int i=0;i<m_seqLenElPtr->m_LenArry.size();i++)
            // {
            //     fwrite(ptr, 1, m_seqLenElPtr->m_LenArry[i], f);
            //     fwrite("\n", 1, 1, f);
            //     ptr += m_seqLenElPtr->m_LenArry[i];
            // }
            // fclose(f);
            arc_stdout(level::ERROR, "seq md5 check fail");
        }
    }
    return true;
}

int SeqProcess::getVal(uint8_t *&ptr, int count)
{
    int ret = 0;
    int tmp;
    for(int i=0;i<count;i++)
    {
        tmp = ptr[i];
        ret += tmp<<i;
    }
    ptr += count;
    return ret;
}

/*-----------------------independent-----------------------*/

int SeqProcess::docompressSeq_independent(char *outptr)
{
    char *ptmp = outptr;

    int szdegetip = compressDegeTip(outptr);
    outptr += szdegetip;
    int szdege = compressDegeCh(outptr);
    outptr += szdege;

    if(m_paramPtr->m_alignType)
    {
        int szorder = compressRefidx(outptr);
        outptr += szorder;
        if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            int szinsr = compressPEInsrNum(outptr);
            outptr += szinsr;
            int szpeRel = compressPERelation(outptr);
            outptr += szpeRel;
        }
        int szpos = compressAlignInfo_Pos(outptr);
        outptr += szpos;
        int szmis = compressAlignInfo_Mis(outptr);
        outptr += szmis;
        int szrev = compressAlignInfo_Rev(outptr);
        outptr += szrev;
        int szcigal = compressAlignInfo_CigaL(outptr);
        outptr += szcigal;
        int szcigav = compressAlignInfo_CigaV(outptr);
        outptr += szcigav;
    }

    int szunmapseq = compressUnmapSeq_independent(outptr);
    outptr += szunmapseq;

    return outptr - ptmp;
}

int SeqProcess::compressUnmapSeq_independent(char *outptr)//没有映射到参考序列
{
    int idlen = Encap::setID(BLOCKSEQ_UNMAPSEQ, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;

    int count = m_seqLenElPtr->m_LenArry.size();    
    RangeCoder rc;
    rc.output(outptr);
    rc.StartEncode();

    for(int i=0;i<m_seq_size;i++)
    {
        m_model_seq[i].reset();
    }

    SIMPLE_MODEL<2> seq_indicate; 
    
    char *seq_p = m_seqElPtr->m_rawBuf.getdata();
    if(m_paramPtr->m_alignType == ALIGNTYPE_REFFREE)
    {
        for (uint32_t i = 0; i < count; i++) 
        {
            
            if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
            {
                encode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
            }
            else
            {
                encode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
            }
            seq_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    else
    {
        uint32_t i_len = 0;
        uint32_t i = 0;
        for (; i < m_seqElPtr->m_vecRefIdx.size(); i++, i_len++) 
        {
            if (m_seqElPtr->m_vecRefIdx[i] == 0)
            {
                uint32_t len = m_seqLenElPtr->m_LenArry[i_len];
                if (len > m_paramPtr->m_MinAlignLen && len < m_paramPtr->m_MaxAlignLen) // 分段
                {
                    uint32_t len1 = len / 2;
                    uint32_t len2 = len - len / 2;
                    char *pseq1 = seq_p;
                    char *pseq2 = seq_p + len1;
                    i++;
                    c_i_divide(pseq1, len1, i, seq_indicate, &rc);
                    i++;
                    c_i_divide(pseq2, len2, i, seq_indicate, &rc);
                }
                else
                {
                    if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
                    {
                        encode_seq_dege(&rc, seq_p, len, seq_indicate);
                    }
                    else
                    {
                        encode_seq_normal(&rc, seq_p, len);
                    }
                }
            }
            // c_i_divide(seq_p, m_seqLenElPtr->m_LenArry[i_len], i, seq_indicate, &rc);
            
            seq_p += m_seqLenElPtr->m_LenArry[i_len];
        }
        
        for(; i_len<count;i_len++, i++)
        {
            if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
            {
                encode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i_len], seq_indicate);
            }
            else
            {
                encode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i_len]);
            }
            seq_p += m_seqLenElPtr->m_LenArry[i_len];
        }

        //--------------------------------------------------
        // uint32_t i = 0;
        // for(; i<m_seqElPtr->m_vecRefIdx.size();i++)
        // {
        //     if(m_seqElPtr->m_vecRefIdx[i] == 0)
        //     {
        //         if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
        //         {
        //             encode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
        //         }
        //         else
        //         {
        //             encode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
        //         }
        //     }
        //     seq_p += m_seqLenElPtr->m_LenArry[i];
        // }
        // for(; i<count;i++)
        // {
        //     if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
        //     {
        //         encode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
        //     }
        //     else
        //     {
        //         encode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
        //     }
        //     seq_p += m_seqLenElPtr->m_LenArry[i];
        // }
        //--------------------------------------------------
    }
    
    rc.FinishEncode();
    int sz = rc.size_out();
    Encap::setSize(sz, SIZENUM4, psize);

    return idlen+SIZENUM4+sz;
}

void SeqProcess::c_i_divide(char *seq_p, uint32_t len, uint32_t &i, SIMPLE_MODEL<2>& seq_indicate, RangeCoder *rc)
{
    if (m_seqElPtr->m_vecRefIdx[i] == 0)
    {
        if (len > m_paramPtr->m_MinAlignLen) // 分段
        {
            uint32_t len1 = len / 2;//不断分割长度进行调用
            uint32_t len2 = len - len / 2;
            char *pseq1 = seq_p;
            char *pseq2 = seq_p + len1;
            i++;
            c_i_divide(pseq1, len1, i, seq_indicate, rc);
            i++;
            c_i_divide(pseq2, len2, i, seq_indicate, rc);
        }
        else
        {
            if(m_seqElPtr->m_vecDegeTip[i]) //存在简并碱基
            {
                encode_seq_dege(rc, seq_p, len, seq_indicate);
            }
            else
            {
                encode_seq_normal(rc, seq_p, len);
            }
        }
    }
}




void SeqProcess::encode_seq_normal(RangeCoder *rc, char *seq, int len) 
{
    int last = 0x007616c7 & NS_MASK;

    uint8_t b = 0;
    for (int i = 0; i < len; i++) 
    {
        b = hash_std_table[seq[i]];
        m_model_seq[last].encodeSymbol(rc, b);
        last = ((last << 2) + b) & NS_MASK; //10个上下文，因为要滑动，*4为了腾出一个空间给新的碱基
    }
}

void SeqProcess::encode_seq_dege(RangeCoder *rc, char *seq, int len, SIMPLE_MODEL<2> &seq_indicate) 
{
    int last = 0x007616c7 & NS_MASK;

    uint8_t b = 0;
    for (int i = 0; i < len; i++) {
        b = hash_std_table[seq[i]];

        if (b > 3) //简并碱基
        {
            seq_indicate.encodeSymbol(rc, 1);
        } else {
            seq_indicate.encodeSymbol(rc, 0);
            m_model_seq[last].encodeSymbol(rc, b);
            last = ((last << 2) + b) & NS_MASK;
        }
    }
}

int SeqProcess::decompressUnmapSeq_independent(char* data)
{
    int count = m_seqLenElPtr->m_LenArry.size();
    RangeCoder rc;
    rc.input(data);
    rc.StartDecode();

    for(int i=0;i<m_seq_size;i++){
        m_model_seq[i].reset();
    }
    
    SIMPLE_MODEL<2> seq_indicate; 

    char *seq_p = m_seqElPtr->m_rawBuf.getdata();
    m_degechidx = 0;
    char *ptmp = seq_p;

    if(m_paramPtr->m_alignType == ALIGNTYPE_REFFREE)
    {
        for (uint32_t i = 0; i < count; i++) {
            if(m_seqElPtr->m_vecDegeTip[i])
            {
                decode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
            }
            else
            {
                decode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
            }
            seq_p += m_seqLenElPtr->m_LenArry[i];
        }
    }
    else
    {
        // uint32_t i = 0, limit;
        // uint64_t realpos = 0;
        // bool ret = false;
        m_posPtr = m_seqElPtr->m_vecBlockPos.getdata();
        m_cigalPtr = m_seqElPtr->m_vecCigaL.getdata();
        m_misptr = m_seqElPtr->m_vecMis.getdata();
        m_mapidx = 0;
        m_cigavidx = 0;
        m_peRelidx = 0;


        uint32_t i_len = 0;
        uint32_t i = 0;
        for (; i < m_seqElPtr->m_vecRefIdx.size(); i++, i_len++) 
        {
            // d_i_divide(seq_p, m_seqLenElPtr->m_LenArry[i_len], i, seq_indicate, &rc);
            //---------
            uint32_t len = m_seqLenElPtr->m_LenArry[i_len];
            if(m_seqElPtr->m_vecRefIdx[i])
            {
                uint32_t limit = len*m_paramPtr->m_iMaxmis/100;
                uint64_t realpos = getRealPos(false, m_seqElPtr->m_vecRefIdx[i]);
                AlignInfoToSeq_independent(m_seqElPtr->m_vecDegeTip[i], seq_p, len, limit, realpos);
            }
            else 
            {
                if (len > m_paramPtr->m_MinAlignLen && len < m_paramPtr->m_MaxAlignLen) // 分段
                { 
                    uint32_t len1 = len / 2;
                    uint32_t len2 = len - len / 2;
                    char *pseq1 = seq_p;
                    char *pseq2 = seq_p + len1;
                    i++;
                    d_i_divide(pseq1, len1, i, seq_indicate, &rc);
                    i++;
                    d_i_divide(pseq2, len2, i, seq_indicate, &rc);
                }
                else
                {
                    if(m_seqElPtr->m_vecDegeTip[i])
                    {
                        decode_seq_dege(&rc, seq_p, len, seq_indicate);
                    }
                    else
                    {
                        decode_seq_normal(&rc, seq_p, len);
                    }
                }
            }
            //---------
            seq_p += m_seqLenElPtr->m_LenArry[i_len];
        }

        for (;i_len < count;i_len++, i++) {
            if(m_seqElPtr->m_vecDegeTip[i])
            {
                decode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i_len], seq_indicate);
            }
            else
            {
                decode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i_len]);
            }
            seq_p += m_seqLenElPtr->m_LenArry[i_len];
        }

        //-------------------------------------------------------
        // for (;i < m_seqElPtr->m_vecRefIdx.size(); i++) {
        //     if(m_seqElPtr->m_vecRefIdx[i])
        //     {
        //         ret = false;
        //         limit = m_seqLenElPtr->m_LenArry[i]*m_paramPtr->m_iMaxmis/100;
        //         if(m_paramPtr->m_arcType == ARCTYPE_PE && (i & 1) && m_seqElPtr->m_vecRefIdx[i-1]>0)
        //         {
        //             ret = true;
        //         }
        //         realpos = getRealPos(ret, m_seqElPtr->m_vecRefIdx[i]);
        //         AlignInfoToSeq_independent(m_seqElPtr->m_vecDegeTip[i], seq_p, m_seqLenElPtr->m_LenArry[i], limit, realpos);
        //     }
        //     else
        //     {
        //         if(m_seqElPtr->m_vecDegeTip[i])
        //         {
        //             decode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
        //         }
        //         else
        //         {
        //             decode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
        //         }
        //     }
            
        //     seq_p += m_seqLenElPtr->m_LenArry[i];
        // }

        // for (;i < count; i++) {
        //     if(m_seqElPtr->m_vecDegeTip[i])
        //     {
        //         decode_seq_dege(&rc, seq_p, m_seqLenElPtr->m_LenArry[i], seq_indicate);
        //     }
        //     else
        //     {
        //         decode_seq_normal(&rc, seq_p, m_seqLenElPtr->m_LenArry[i]);
        //     }
        //     seq_p += m_seqLenElPtr->m_LenArry[i];
        // }
        //-------------------------------------------------------
    }
    
    rc.FinishDecode();
    int size = seq_p - ptmp;
    m_seqElPtr->m_rawBuf.updatesize(size);
    return size; 
}

// 组装出seq，递归
void SeqProcess::d_i_divide(char *seq_p, uint32_t len, uint32_t &i, SIMPLE_MODEL<2>& seq_indicate, RangeCoder *rc)
{
        // uint32_t i = 0, limit;
        // uint64_t realpos = 0;
        // bool ret = false;
    if(m_seqElPtr->m_vecRefIdx[i])
    {
        uint32_t limit = len*m_paramPtr->m_iMaxmis/100;
        uint64_t realpos = getRealPos(false, m_seqElPtr->m_vecRefIdx[i]);
        AlignInfoToSeq_independent(m_seqElPtr->m_vecDegeTip[i], seq_p, len, limit, realpos);
    }
    else 
    {
        if (len > m_paramPtr->m_MinAlignLen) // 分段
        { 
            uint32_t len1 = len / 2;
            uint32_t len2 = len - len / 2;
            char *pseq1 = seq_p;
            char *pseq2 = seq_p + len1;
            i++;
            d_i_divide(pseq1, len1, i, seq_indicate, rc);
            i++;
            d_i_divide(pseq2, len2, i, seq_indicate, rc);
        }
        else
        {
            if(m_seqElPtr->m_vecDegeTip[i])
            {
                decode_seq_dege(rc, seq_p, len, seq_indicate);
            }
            else
            {
                decode_seq_normal(rc, seq_p, len);
            }
        }
    }
}

void SeqProcess::decode_seq_normal(RangeCoder *rc, char *seq, int len) 
{
    int last = 0x007616c7 & NS_MASK;

    uint8_t b = 0;
    for (int i = 0; i < len; i++) {
        b = m_model_seq[last].decodeSymbol(rc);
        last = ((last << 2) + b) & NS_MASK;
        seq[i] = fqzdec[b];
    }
}

void SeqProcess::decode_seq_dege(RangeCoder *rc, char *seq, int len, SIMPLE_MODEL<2> &seq_indicate) 
{
    int last = 0x007616c7 & NS_MASK;

    uint8_t b = 0;
    for (int i = 0; i < len; i++) {
        if (seq_indicate.decodeSymbol(rc)) //简并碱基
        {
            seq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
        } else {
            b = m_model_seq[last].decodeSymbol(rc);
            last = (last * 4 + b) & NS_MASK;
            seq[i] = fqzdec[b];
        }
    }
}

void SeqProcess::AlignInfoToSeq_independent(bool bdege, char *pseq, int seqlen, int limit, uint64_t realpos)
{
    bool isrev = m_seqElPtr->m_vecRev[m_mapidx++];
    // int cigar_num = getVal(m_misptr,getbitnum(limit));
    // decodeAlignInfo(seqlen, cigar_num, realpos);
    decodeAlignInfo(seqlen, 0, realpos);

    if(bdege)
    {
        if(isrev)
        {
            int j = seqlen - 1;
            for (int i = 0; i < seqlen; i++, j--) {
                if(m_refseq[j] == 4)
                {
                    pseq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                }
                else
                {
                    pseq[i] = "ACGT"[m_refseq[j] ^ 3];
                }
            }
        }
        else
        {
            for (int i = 0; i < seqlen; i++) {
                if(m_refseq[i] == 4)
                {
                    pseq[i] = m_seqElPtr->m_vecDegeCh[m_degechidx++];
                }
                else
                {
                    pseq[i] = "ACGT"[m_refseq[i]&3];
                }
                
            }
        }
    }
    else
    {
        if(isrev)
        {
            int j = seqlen - 1;
            for (int i = 0; i < seqlen; i++, j--) {
                pseq[i] = "ACGT"[m_refseq[j] ^ 3];
            }
        }
        else
        {
            for (int i = 0; i < seqlen; i++) {
                pseq[i] = "ACGT"[m_refseq[i]&3];
            }
        }
    }
}

