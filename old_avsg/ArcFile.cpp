#include "ArcFile.h"
#include <fcntl.h> 

ArcFile::ArcFile()
{
    m_paramPtr = Param::GetInstance();
    if(m_paramPtr->m_actionType == ACTIONTYPE_DOENCODE)
    {
        m_paramPtr->m_outfd[0] = open(m_paramPtr->m_strArcFile.c_str(), O_RDWR|O_CREAT, 0664);
        if(m_paramPtr->m_outfd[0] < 0)
        {
            arc_stdout(level::ERROR, "arc file create fail");
        }
    }
    if(m_paramPtr->m_actionType == ACTIONTYPE_DODECODE || m_paramPtr->m_actionType == ACTIONTYPE_EXTRACT )
    {
        m_filesize = getFileSize(m_paramPtr->m_strArcFile.c_str());
        m_file.open(m_paramPtr->m_strArcFile, std::ios::binary | std::ios::in);
        m_file.read(m_readbuf, 64);
        readInfoEncap();        
    }
}

ArcFile::~ArcFile()
{
    m_file.close();
}

int ArcFile::writeParam()
{
    char *ptmp = m_bufPtr;
    int idlen = Encap::setID(ARCTAIL_PARAM, m_bufPtr);
    m_bufPtr += idlen;
    char *pbuf = m_bufPtr;
    m_bufPtr += 4;
    int sz = m_paramPtr->ParamEncap(m_bufPtr);
    Encap::setSize(sz, 4, pbuf);
    m_bufPtr += sz;
    return m_bufPtr - ptmp;
}

uint64_t ArcFile::writeFileInfo()
{
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        m_paramPtr->m_iBlocknum = m_paramPtr->m_vecBlockInfoPE.size();
        for (uint32_t i = 0; i < m_paramPtr->m_iBlocknum; i++) {
            m_totalBlockSize += m_paramPtr->m_vecBlockInfoPE[i].compressSize;
        }
    }
    else
    {
        m_paramPtr->m_iBlocknum = m_paramPtr->m_vecBlockInfoSE.size();
        for (uint32_t i = 0; i < m_paramPtr->m_iBlocknum; i++) {
            m_totalBlockSize += m_paramPtr->m_vecBlockInfoSE[i].compressSize;
        }        
    }

    //write head
    char headbuf[32]={0};
    char *pbuf = headbuf;
    memcpy(pbuf, ".arc", 4);
    pbuf += 4;

    *pbuf++ = (MAJOR_VERS) & 0xff;
    *pbuf++ = (MINOR_VERS) & 0xff;
    *pbuf++ = (STAGE_VERS) & 0xff;
    
    int len = Encap::setID(ARCPART_BLOCK, pbuf);
    pbuf += len;
    Encap::setSize(m_totalBlockSize, 8, pbuf);
    pbuf += 8;
    int headlen = pbuf - headbuf;
    if(headlen != ARCHEADLEN)
    {
        arc_stdout(level::ERROR, "headlen err %d", headlen);
    }
    pwrite(m_paramPtr->m_outfd[0], headbuf, headlen, 0);

    //write tail
    m_bufPtr = m_tailbuf.mem_malloc(m_paramPtr->m_iBlockSize);
    int idlen = Encap::setID(ARCPART_TAIL, m_bufPtr);
    m_bufPtr += idlen+4;

    int len1 = writeParam();

    int taillen = m_bufPtr - m_tailbuf.getdata();
    char *psize = m_tailbuf.getdata() + idlen;
    int size = taillen-idlen-4;
    Encap::setSize(size, 4, psize);

    pwrite(m_paramPtr->m_outfd[0], m_tailbuf.getdata(), taillen, headlen+m_totalBlockSize);
    m_filesize = headlen+m_totalBlockSize+taillen;
    if(close(m_paramPtr->m_outfd[0]) != 0)
    {
        arc_stdout(level::ERROR,"close error");
    }

    return m_filesize;
}

int ArcFile::readInfoEncap()
{
    char *pbuf = m_readbuf;
    if(memcmp(pbuf, ".arc", 4)==0)
    {
        pbuf += 4;
        int major = *pbuf++;
        int minor = *pbuf++;
        int stage = *pbuf++;
        EncapInfo info;
        Encap::getValue(pbuf, info);
        if(info.val != ARCPART_BLOCK)
        {
            arc_stdout(level::ERROR, "arc file is old version");
        }
        pbuf += info.len;
        Encap::getValue(pbuf, info);
        m_tailoffset += ARCHEADLEN;
        m_tailoffset += info.val;
    }
    else
    {
        arc_stdout(level::ERROR, "File is not arc type");
    }

    return 0;
}

int ArcFile::readBlockParam(std::map<uint32_t, BlockParamInfo> &mapBParam)
{
    char *pbuf = m_readbuf;
    EncapInfo info;
    uint64_t  offset = ARCHEADLEN;
    pbuf += ARCHEADLEN;
    BlockParamInfo param;
    char buf[64];
    int i = 0;
    while (offset < m_tailoffset)
    {
        pbuf += 1;
        Encap::getValue(pbuf, info);
        offset += info.val + 5;

        pbuf += 11;
        param.onech = (uint8_t)pbuf[0];
        pbuf += 3;
        param.blockidx = DECODE_INT((unsigned char*) pbuf);
        pbuf += 6;
        param.readcount = DECODE_INT((unsigned char*) pbuf);
        param.writeidx = i++;

        mapBParam.emplace(param.blockidx, param);

        m_file.seekg(offset);
        m_file.read(buf, 64);
        pbuf = buf;
    }
    return 0;
}

int ArcFile::readTail()
{
    EncapInfo info;
    int taillen = m_filesize - m_tailoffset;
    m_bufPtr = m_tailbuf.mem_malloc(taillen);
    m_file.seekg(m_tailoffset);
    m_file.read(m_bufPtr, taillen);

    Encap::getValue(m_bufPtr, info);
    if(info.val != ARCPART_TAIL)
    {
        arc_stdout(level::ERROR, "readtail type err");
    }

    m_bufPtr += info.len;
    Encap::getValue(m_bufPtr, info);
    taillen = info.val;
    m_bufPtr += info.len;
    
    uint32_t r_len = 0;
    while (r_len < taillen)
    {
        Encap::getValue(m_bufPtr, info);
        int id = info.val;
        m_bufPtr += info.len;
        r_len += info.len;
        Encap::getValue(m_bufPtr, info);
        m_bufPtr += info.len;
        r_len += info.len;

        doReadTail(id, info.val, m_bufPtr);
        m_bufPtr += info.val;
        r_len += info.val;
    }
    return 0;
}


int ArcFile::doReadTail(int id, uint32_t size, char* data)
{
    switch (id)
    {
    case ARCTAIL_PARAM:
        m_paramPtr->ParamUnEncap(size, data);
        break;
    default:
        break;
    }

    return 0;
}
