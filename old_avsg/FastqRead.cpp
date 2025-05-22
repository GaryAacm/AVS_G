#include "FastqRead.h"
#include "BufPool.h"
#include "Param.h"
#include "IRef.h"
#include "ThreadPool.h"
#include <zlib.h>
#include "EncodeSEWorker.h"
#include "EncodePEWorker.h"
#include "log/progress_bar.h"


FastqRead::FastqRead(IRef *refPtr):
    m_refPtr(refPtr)
{
    m_paramPtr = Param::GetInstance();
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        m_vecpos1.mem_malloc(1024*1024);
        m_vecpos2.mem_malloc(1024*1024);
    }
    m_memBufPoolPtr = m_paramPtr->m_memBufPoolPtr;
    m_thPool = new ThreadPool();
}

FastqRead::~FastqRead()
{
    RELEASEPTR(m_thPool);
}

bool FastqRead::isFinish()
{
    return m_thPool->isAllFinish();
}

int FastqRead::getProcessCount()
{
    return m_thPool->getProcessCount();
}

int FastqRead::getEndPos(int start, char *pbuf)
{
    int samenum = 0;
    for(int i=start;i>0;i--)
    {
        if(pbuf[i] == '\n' && pbuf[i+1] == '@') //第一行或者第四行
        {
            for(int j=0, k= i+1;j < m_strfirstline.length();j++, k++)
            {
                if(m_strfirstline[j] == pbuf[k])
                {
                    samenum++;
                    if(samenum > 5) //判定是第一行
                    {
                        return i;
                    }
                }
                else
                {
                    samenum = 0;
                }
            }

        }
    }
    arc_stdout(level::ERROR, "not find endpos");
    return 0;
}

void FastqRead::cultPEbuf(MemBuf *ptr, bool bread)
{
    char *pbuf1 = ptr->m_bufptr[0];
    char *pbuf2 = ptr->m_bufptr[1];
    int readlen1 = ptr->m_len[0];
    int readlen2 = ptr->m_len[1];
    m_vecpos1.clear();
    m_vecpos2.clear();
    int i=0,j=0;
    for(;i<readlen1 && j< readlen2;i++,j++)
    {
        if(pbuf1[i] == '\n')
        {
            m_vecpos1.push_back(i);
        }
        if(pbuf2[j] == '\n')
        {
            m_vecpos2.push_back(j);
        }
    }

    while (i<readlen1)
    {
        if(pbuf1[i] == '\n')
        {
            m_vecpos1.push_back(i);
        }
        i++;
    }

    while (j< readlen2)
    {
        if(pbuf2[j] == '\n')
        {
            m_vecpos2.push_back(j);
        }
        j++;
    }

    int min = m_vecpos1.size() < m_vecpos2.size() ? m_vecpos1.size() : m_vecpos2.size();
    if(bread)
    {
        min -= 2;
        int start = m_vecpos1[min];

        int endpos = getEndPos(start, pbuf1);
        while (m_vecpos1[min] != endpos)
        {
            min--;
        }
        if(endpos != m_vecpos1[min])
        {
            arc_stdout(level::ERROR, "pe cult err");
        }
        ptr->m_len[0] = m_vecpos1[min] +1;
        ptr->m_len[1] = m_vecpos2[min] +1;
        m_leftstr[0].append(&pbuf1[ptr->m_len[0]], readlen1-ptr->m_len[0]);
        m_leftstr[1].append(&pbuf2[ptr->m_len[1]], readlen2-ptr->m_len[1]);
    }
    else
    {
        ptr->m_len[0] = readlen1;
        ptr->m_len[1] = readlen2;
    }

    ptr->m_vecpos1.swap(m_vecpos1);
    ptr->m_vecpos2.swap(m_vecpos2);
    ptr->m_penum = min;
}


void FastqRead::cultbuf(MemBuf *ptr)
{
    char *pbuf = ptr->m_bufptr[0];
    int start = ptr->m_len[0]-m_strfirstline.length();
    int readlen = ptr->m_len[0];
    int endpos = getEndPos(start, pbuf);
    ptr->m_len[0] = endpos+1;
    m_leftstr[0].append(&pbuf[ptr->m_len[0]], readlen-ptr->m_len[0]);
}

void FastqRead::startEncodeThread()
{
    int thnum = m_thPool->getThreadNum();
    if(thnum < m_memBufPoolPtr->getFullbufNum()&&
       thnum < m_paramPtr->m_iThreadnum)
    {
        IWorker *ptr = nullptr;
        if(m_paramPtr->m_arcType)
        {
            ptr = new EncodeSEWorker(thnum, m_refPtr);
        }
        else
        {
            ptr = new EncodePEWorker(thnum, m_refPtr);
        }
        ptr->init();//初始化获取数据和进程
        m_thPool->addJob(ptr);
    }
}

void FastqRead::sendEndMsg()
{
    int thnum = m_thPool->getThreadNum();
    for(int i=0;i<thnum;i++)
    {
        MemBuf *ptr = m_memBufPoolPtr->getEmptybuf(false);
        ptr->m_len[0] = 0;
        m_memBufPoolPtr->addFullbuf(ptr);
    }
    //arc_stdout(level::INFO, "threadnum:%d bufnum:%d", thnum, m_memBufPoolPtr->getTotalbufNum());
}

void FastqRead::getFirstLine(char *pbuf)
{
    if(m_strfirstline.empty())
    {
        int pos[4]={0};
        int i=0,j=0;
        while (i<4)
        {
            if(pbuf[j++] == '\n')
            {
                pos[i++] = j;
            }
        }

        m_strfirstline.append(pbuf, pos[0]);
        if(pos[2]-pos[1] > 2)
        {
            m_bonech = false;
        }
    }
}

void FastqRead::ReadFastq()
{
    auto rtime = Profile::getInstance().getModule("main")->getTimer();
    rtime->start("ReadTime");
    uint64_t offset = 0;
    int readlen=0, leftlen = 0, reallen=0;
    bool bread = true;
    int fileidx = 0;
    for(std::string &strFile : m_paramPtr->m_vecFqFile)
    {
        gzFile file = gzopen(strFile.c_str(), "r");
        gzbuffer(file, 256*1024);
        readlen = 0, leftlen = 0, reallen = 0;
        m_strfirstline.clear();
        bread = true;
        m_bonech = true;
        offset = 0;
        int blockidx = 0;
        while (bread)//为压缩的准备
        {
            MemBuf *ptr = m_memBufPoolPtr->getEmptybuf();
            char *pbuf = ptr->m_bufptr[0];
            leftlen = m_leftstr[0].length();

            memcpy(pbuf, m_leftstr[0].c_str(), leftlen);

            m_leftstr[0].clear();
            pbuf += leftlen;
            readlen = m_paramPtr->m_iBlockSize-leftlen;
            reallen = gzread(file, pbuf, readlen);//存储在pbuf里面

            if(reallen == -1) //error
            {
                int z_errnum = 0;
                const char *errmsg = gzerror(file, &z_errnum);
                arc_stdout(level::ERROR, "se read error %s", errmsg);
            }

            getFirstLine(pbuf);

            ptr->m_len[0] = reallen + leftlen;
            if(reallen < readlen) //end
            {
                bread = false;
            }
            else
            {
                cultbuf(ptr);
            }

            ptr->m_offset[0] = offset;
            ptr->m_fileidx = fileidx;
            ptr->m_bonech = m_bonech;
            ptr->m_blockidx = blockidx;
            blockidx++;
            offset += ptr->m_len[0];
            m_memBufPoolPtr->addFullbuf(ptr);

            startEncodeThread();
            m_blocknum++;
            m_barPtr->update(m_blocknum/3);
        }

        gzclose(file);
        fileidx++;
    }
    rtime->stop();
    sendEndMsg();
}


void FastqRead::ReadPEFastq()
{
    auto rtime = Profile::getInstance().getModule("main")->getTimer();
    rtime->start("ReadTime");
    uint64_t offset[2] = {0};
    gzFile file[2];
    file[0] = gzopen(m_paramPtr->m_vecFqFile[0].c_str(), "r");
    file[1] = gzopen(m_paramPtr->m_vecFqFile[1].c_str(), "r");
    gzbuffer(file[0], 256*1024);
    gzbuffer(file[1], 256*1024);
    int readlen[2]={0}, reallen[2]={0};
    int leftlen = 0;
    m_strfirstline.clear();
    char *pbuf = nullptr;
    int halfsize = m_paramPtr->m_iBlockSize/2;
    bool bread = true;
    int blockidx = 0;
    while (bread)
    {

        MemBuf *ptr = m_memBufPoolPtr->getEmptybuf();

        for(int i=0;i<2;i++)
        {
            pbuf = ptr->m_bufptr[i];
            leftlen = m_leftstr[i].length();
            memcpy(pbuf, m_leftstr[i].c_str(), leftlen);
            m_leftstr[i].clear();
            pbuf += leftlen;
            readlen[i] = halfsize-leftlen;
            reallen[i] = gzread(file[i], pbuf, readlen[i]);
            if(reallen[i] == -1) //error
            {
                int z_errnum = 0;
                const char *errmsg = gzerror(file[i], &z_errnum);
                arc_stdout(level::ERROR, "pe read error %d %s", i, errmsg);
            }
            ptr->m_len[i] = reallen[i] + leftlen;
        }

        getFirstLine(pbuf);

        if(reallen[0] < readlen[0] &&
           reallen[1] < readlen[1] )
        {
            bread = false;
        }

        cultPEbuf(ptr, bread);

        ptr->m_bonech = m_bonech;
        ptr->m_offset[0] = offset[0];
        offset[0] += ptr->m_len[0];
        ptr->m_offset[1] = offset[1];
        ptr->m_blockidx = blockidx;
        blockidx++;
        offset[1] += ptr->m_len[1];
        m_memBufPoolPtr->addFullbuf(ptr);

        startEncodeThread();
        m_blocknum++;
        m_barPtr->update(m_blocknum/3);
    }

    gzclose(file[0]);
    gzclose(file[1]);
    rtime->stop();
    sendEndMsg();
}


void FastqRead::ReadPipe()
{
    uint64_t offset = 0;
    int readlen=0, leftlen = 0, reallen=0;
    m_strfirstline.clear();
    bool bread = true;
    int blockidx = 0;
    while (bread)
    {
        MemBuf *ptr = m_memBufPoolPtr->getEmptybuf();

        char *pbuf = ptr->m_bufptr[0];
        leftlen = m_leftstr[0].length();
        memcpy(pbuf, m_leftstr[0].c_str(), leftlen);
        m_leftstr[0].clear();
        pbuf += leftlen;
        readlen = m_paramPtr->m_iBlockSize-leftlen;
        reallen = fread(pbuf, 1, readlen, stdin);
        if(reallen == -1) //error
        {
            arc_stdout(level::ERROR, "stdin error");
        }

        getFirstLine(pbuf);

        ptr->m_len[0] = reallen + leftlen;
        if(reallen < readlen)
        {
            arc_stdout(level::ERROR, "stdin end");
            bread = false;
        }
        else
        {
            cultbuf(ptr);
        }

        ptr->m_offset[0] = offset;
        ptr->m_bonech = m_bonech;
        ptr->m_blockidx = blockidx;
        blockidx++;
        offset += ptr->m_len[0];
        m_memBufPoolPtr->addFullbuf(ptr);

        startEncodeThread();
    }

    sendEndMsg();
}