#include "ThreadPool.h"

std::mutex IWorker::m_mtx_BlockInfo;
ThreadPool::ThreadPool()
{
}

ThreadPool::~ThreadPool()
{
    for(int i=0;i<m_vecThread.size();i++)
    {
        Thread *thPtr = m_vecThread[i];
        RELEASEPTR(thPtr);
    }
}

void ThreadPool::addJob(IWorker *workPtr)
{
    assert(workPtr);
    Thread *thPtr = new Thread(workPtr);
    assert(thPtr);
    m_vecThread.emplace_back(thPtr);
    thPtr->startThread();
    m_threadnum++;
}

int ThreadPool::getThreadNum()
{
    return m_threadnum;
}

bool ThreadPool::isAllFinish()
{
    if(m_vecThread.empty()) return false;
    for(int i=0;i<m_vecThread.size();i++)
    {
        Thread *thPtr = m_vecThread[i];
        if(thPtr->getstate() != THREAD_END)
        {
            return false;
        }
    }
    return true;
}

int ThreadPool::getProcessCount()
{
    int count = 0;
    for(int i=0;i<m_vecThread.size();i++)
    {
        Thread *thPtr = m_vecThread[i];
        count += thPtr->getProcessCount();
    }
    return count;
}