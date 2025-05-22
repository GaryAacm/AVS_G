#include "Param.h"

Param::Param()
{
}

Param::~Param()
{
    RELEASEPTR(m_memBufPoolPtr);
    close(m_outfd[0]);
    close(m_outfd[1]);
}

void Param::add(const string &str)
{
    m_vecAddParam.emplace_back(str);
}

int Param::EncapRefMd5(uint32_t id, char *pbuf)
{
    if(m_alignType)
    {
        getRefMd5();
        char *ptmp = pbuf;
        int idlen = Encap::setID(id, pbuf);
        pbuf += idlen;
        Encap::setSize(16, 1, pbuf);
        pbuf += 1;
        memcpy(pbuf, m_refmd5, 16);
        pbuf += 16;
        return pbuf - ptmp;  
    }
    return 0;
}


void Param::getRefMd5()
{
    fstream md5_s;
    char md5_path[PATH_LEN] = {0};
    sprintf(md5_path, "%s.md5", m_strRefFile.c_str());
    md5_s.open(md5_path, std::ios::binary | std::ios::in);
    md5_s.read((char*)m_refmd5, 16);
    md5_s.close();
}

int Param::EncapBlockInfoSE(char *pbuf)
{
    int len = sizeof(BlockInfoSE);
    int size = len*m_iBlocknum;
    
    for (uint32_t i = 0; i < m_iBlocknum; i++) {
        memcpy(pbuf, &m_vecBlockInfoSE[i], len);
        pbuf += len;
    }
    return size;
}

int Param::EncapBlockInfoPE(char *pbuf)
{
    int len = sizeof(BlockInfoPE);
    int size = len*m_iBlocknum;

    for (uint32_t i = 0; i < m_iBlocknum; i++) {
        memcpy(pbuf, &m_vecBlockInfoPE[i], len);
        pbuf += len;
    }
    
    return size;
}

int Param::EncapBlockInfo(uint32_t id, char *pbuf)
{
    int idlen = Encap::setID(id, pbuf);
    pbuf += idlen;
    char *psize = pbuf;
    pbuf += 4;
    int sz = 0;
    if(m_arcType == ARCTYPE_PE)
    {
        sz = EncapBlockInfoPE(pbuf);
    }
    else
    {  
        sz = EncapBlockInfoSE(pbuf);
    }   

    Encap::setSize(sz, 4, psize); 
    return sz + 4 + idlen;
}

int Param::ParamEncap(char *bufPtr)
{
    char *pbuf = bufPtr;
    int sz = Encap::encapUint8(PARAM_DATACHECK, m_checktype, pbuf);
    pbuf += sz;
    sz = Encap::encapUint8(PARAM_ALIGNTYPE, m_alignType, pbuf);
    pbuf += sz;
    sz = Encap::encapUint8(PARAM_ARCTYPE, m_arcType, pbuf);
    pbuf += sz;
    sz = Encap::encapUint8(PARAM_DEPENDENT, m_dependent, pbuf);
    pbuf += sz;
    sz = Encap::encapUint32(PARAM_IBLOCKNUM, m_iBlocknum, pbuf);
    pbuf += sz;
    sz = Encap::encapUint32(PARAM_IMAXMIS, m_iMaxmis, pbuf);
    pbuf += sz;
    sz = Encap::encapUint32(PARAM_IBLOCKSIZE, m_iBlockSize, pbuf);
    pbuf += sz;
    sz = Encap::encapUint32(PARAM_KMER_SIZE, m_kmer_size, pbuf);
    pbuf += sz;
    sz = Encap::encapUint32(PARAM_SEQ_COMPRESS_MODE, m_seq_compress_mode, pbuf);
    pbuf += sz;
    sz = EncapRefMd5(PARAM_REFMD5, pbuf);
    pbuf += sz;

    for(string &str : m_vecFqFile)
    {
        getFileName(str);
        sz = Encap::encapStr(PARAM_FQFILE, str, pbuf);
        pbuf += sz;
    }

    for(string &str : m_vecAddParam)
    {
        sz = Encap::encapStr(PARAM_ADDPARAM, str, pbuf);
        pbuf += sz;
    }

    sz = EncapBlockInfo(PARAM_BLOCKINFO, pbuf);
    pbuf += sz;

    sz = pbuf - bufPtr;
    return sz;
}

void Param::ParamUnEncap(uint32_t size, char *buf)
{
    char *pbuf = buf;
    EncapInfo info;
    Encap encap;
    uint32_t rlen = 0;
    while (rlen < size)
    {
        Encap::getValue(pbuf, info);
        int id = info.val;
        pbuf += info.len;
        rlen += info.len;
        Encap::getValue(pbuf, info);
        pbuf += info.len;
        rlen += info.len;

        decompressParam(id, info.val, pbuf); 

        pbuf += info.val;
        rlen += info.val;
    }
}

void Param::decompressParam(int id, int size, char *pbuf)
{
    switch (id)
    {
    case PARAM_DATACHECK:
        m_checktype = (DataCheck)((uint8_t)pbuf[0]);
        break;
    case PARAM_ALIGNTYPE:
        m_alignType = (AlignType)((uint8_t)pbuf[0]);
        break;
    case PARAM_ARCTYPE:
        m_arcType = (ArcType)((uint8_t)pbuf[0]);
        break;
    case PARAM_DEPENDENT:
        m_dependent = (Dependent)((uint8_t)pbuf[0]);
        break;
    case PARAM_IBLOCKNUM:
        m_iBlocknum = DECODE_INT((unsigned char*)pbuf);
        if(m_iThreadnum > m_iBlocknum)
        {
            m_iThreadnum = m_iBlocknum;
        }
        break;
    case PARAM_IMAXMIS:
        m_iMaxmis = DECODE_INT((unsigned char*)pbuf);
        break;
    case PARAM_IBLOCKSIZE:
        m_iBlockSize = DECODE_INT((unsigned char*)pbuf);
        break;
    case PARAM_KMER_SIZE:
        m_kmer_size = DECODE_INT((unsigned char*)pbuf);
        break;
    case PARAM_SEQ_COMPRESS_MODE:
        m_seq_compress_mode = (SeqCompressMode)DECODE_INT((unsigned char*)pbuf);
        break;
    case PARAM_REFMD5:
        {
            getRefMd5();
            if(memcmp(pbuf, m_refmd5, 16)!=0)
            {
                arc_stdout(level::ERROR, "Ref Md5 is not equal");
            }
        }
        break;
    case PARAM_FQFILE:
    {
        string str(pbuf, size);
        m_vecFqFile.emplace_back(str);
    }
        break;
    case PARAM_ADDPARAM:
    {
        string str(pbuf, size);
        m_vecAddParam.emplace_back(str);
    }
        break;
    case PARAM_BLOCKINFO:
        UnEncapBlockInfo(pbuf);
        break;
    default:
        break;
    }
}

void Param::addSEinfo(BlockInfoSE &info)
{
    BlockMetaInfo metaInfo;
    metaInfo.fileidx = info.fileidx;
    metaInfo.compressSize = info.compressSize;
    metaInfo.originalSize1 = info.originalSize;
    metaInfo.originalSize2 = 0;
    metaInfo.compressOffset = info.compressoffset;
    metaInfo.originalOffset1 = info.originaloffset;
    metaInfo.originalOffset2 = 0;

    m_vecBlockMeta.emplace_back(metaInfo);
}

int Param::UnEncapBlockInfoSE(char *pbuf)
{
    int len = sizeof(BlockInfoSE);
    char *ptr = pbuf;
    if(m_arcType == ARCTYPE_SE)
    {
        for(uint32_t i=0;i<m_iBlocknum;i++)
        {
            BlockInfoSE info;
            memcpy(&info, ptr, len);
            ptr += len;

            addSEinfo(info);
        }
    }
    else if(m_arcType == ARCTYPE_MULTI)
    {
        if(!m_hashMultiple.empty())
        {
            for(uint32_t i=0;i<m_iBlocknum;i++)
            {
                BlockInfoSE info;
                memcpy(&info, ptr, len);
                ptr += len;

                if(m_hashMultiple.find(info.fileidx) != m_hashMultiple.end())
                {
                    addSEinfo(info);
                }
            }
        }
        else  //全部解压
        {
            for(uint32_t i=0;i<m_iBlocknum;i++)
            {
                BlockInfoSE info;
                memcpy(&info, ptr, len);
                ptr += len;
                addSEinfo(info);
                m_hashMultiple.insert(info.fileidx);
            }
        }
    }

    return 0;
}

int Param::UnEncapBlockInfoPE(char *pbuf)
{
    int len = sizeof(BlockInfoPE);
    char *ptr = pbuf;
    for(uint32_t i=0;i<m_iBlocknum;i++)
    {
        BlockInfoPE info;
        memcpy(&info, ptr, len);
        ptr += len;

        BlockMetaInfo metaInfo;
        metaInfo.fileidx = 0;
        metaInfo.compressSize = info.compressSize;
        metaInfo.originalSize1 = info.originalSize[0];
        metaInfo.originalSize2 = info.originalSize[1];
        metaInfo.compressOffset = info.compressoffset;
        metaInfo.originalOffset1 = info.originaloffset[0];
        metaInfo.originalOffset2 = info.originaloffset[1];

        m_vecBlockMeta.emplace_back(metaInfo);
    }
    return 0;
}

int Param::UnEncapBlockInfo(char *pbuf)
{
    if(m_arcType == ARCTYPE_PE)
    {
        return UnEncapBlockInfoPE(pbuf);
    }
    else
    {
        return UnEncapBlockInfoSE(pbuf);
    }
}