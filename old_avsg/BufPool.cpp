#include "BufPool.h"
#include "util.h"
#include "Param.h"

MemBuf::MemBuf(bool binit)
{
    if(binit)
    {
        int type = Param::GetInstance()->m_arcType;
        int blocksize = Param::GetInstance()->m_iBlockSize;
        if(Param::GetInstance()->m_PipeOut)
        {
            m_bufptr[0] = new char[blocksize];
        }
        else
        {
            switch (type)
            {
            case ARCTYPE_PE:
            {
                int halfsize = blocksize/2;
                m_bufptr[0] = new char[halfsize];
                m_bufptr[1] = new char[halfsize];
                m_vecpos1.mem_malloc(1024*1024);
                m_vecpos2.mem_malloc(1024*1024); 
            }
                break;
            case ARCTYPE_SE:
            case ARCTYPE_MULTI:
                m_bufptr[0] = new char[blocksize];
                break;
            default:
                break;
            }
        }
    }
}

MemBuf::~MemBuf()
{
    RELEASEARRYPTR(m_bufptr[0]);
    RELEASEARRYPTR(m_bufptr[1]);
}


MemBufPool::MemBufPool(int limit):
    m_limit(limit)
{
}

MemBufPool::~MemBufPool()
{
    for(MemBuf *ptr : m_fullList)
    {
        RELEASEPTR(ptr);
    }
    for(MemBuf *ptr : m_emptyList)
    {
        RELEASEPTR(ptr);
    }
}

MemBuf *MemBufPool::getEmptybuf(bool binit)
{
    MemBuf *ptr = nullptr;
    std::unique_lock<std::mutex> tlock(m_mtx_empty);
    if(m_emptyList.empty())
    {
        if(m_totalbuf < m_limit)
        {

            ptr = new MemBuf(binit);
            if(binit) m_totalbuf++;
            //printf("new buf %d\n", m_totalbuf);
        }
        else //超过了限定,需要等待
        {
            while (m_emptyList.empty())
            {
                //printf("getEmptyJob wait\n");
                m_cv_empty.wait(tlock);
            }

            //printf("%d getEmptyJob wait a buf\n", getpid());
            ptr = m_emptyList.front();
            m_emptyList.pop_front();
            
        }
    }
    else
    {
        //printf("%d getEmptyJob front\n", getpid());
        ptr = m_emptyList.front();
        m_emptyList.pop_front();
    }

    if(ptr == nullptr) 
    {
        ptr = new MemBuf(binit);
        m_totalbuf++;
        printf("add readbuf %d\n", m_totalbuf);
    }

    return ptr;
}

MemBuf *MemBufPool::getFullbuf()
{
    MemBuf *ptr = nullptr;
    std::unique_lock<std::mutex> tlock(m_mtx_full);

    while (m_fullList.empty())
    {
        //printf("%d getFullJob wait\n", num);
        m_cv_full.wait(tlock);
    }

    //printf("%d getFullJob front\n", getpid());
    ptr = m_fullList.front();
    m_fullList.pop_front();
    return ptr;
}

MemBuf *MemBufPool::getIdxFullbuf(int idx)
{
    MemBuf *ptr = nullptr;
    std::unique_lock<std::mutex> tlock(m_mtx_full);

    while (m_fullList.empty())
    {
        m_cv_full.wait(tlock);
    }

    for(auto itor = m_fullList.begin(); itor!=m_fullList.end();itor++)
    {
        if ((*itor)->m_blockidx == idx)
        {
            ptr = *itor;
            m_fullList.erase(itor);
            return ptr;
        }
    }
    return ptr;
}

void MemBufPool::addEmptybuf(MemBuf *ptr)
{
    std::lock_guard<std::mutex> tlock(m_mtx_empty);
    m_emptyList.push_back(ptr);
    m_cv_empty.notify_one();
}

void MemBufPool::addFullbuf(MemBuf *ptr)
{
    std::lock_guard<std::mutex> tlock(m_mtx_full);
    m_fullList.push_back(ptr);
    m_cv_full.notify_one();
}

int MemBufPool::getFullbufNum()
{
    std::lock_guard<std::mutex> tlock(m_mtx_full);
    //printf("fullnum %d\n", m_fullList.size());
    return m_fullList.size();
}




