#ifndef AVS_API_THREADPOOL_H
#define AVS_API_THREADPOOL_H

#include "Param.h"
#include "IWorker.h"
#include <thread>


enum threadstate
{
    THREAD_CREATE = 1,
    THREAD_RUNING,
    THREAD_END
};

class Thread
{
public:
    Thread(IWorker *ptr):m_wroker(ptr){}
    ~Thread(){}
    void startThread()
    {
        m_state = THREAD_CREATE;
        m_thread = std::thread(&Thread::Run, this);
        m_thread.detach();
    }
    void Run()
    {
        m_state = THREAD_RUNING;
        m_wroker->doJob();
        m_state = THREAD_END;
    }
    int getstate()
    {
        return m_state;
    }
    int getProcessCount()
    {
        return m_wroker->getProcessCount();
    }
private:
    int m_state = 0;
    std::thread m_thread;
    IWorker *m_wroker = nullptr;
};

class ThreadPool
{
public:
    ThreadPool();
    ~ThreadPool();
    void addJob(IWorker *workPtr);
    int getThreadNum();
    bool isAllFinish();
    int getProcessCount();
private:
    int m_threadnum = 0;
    std::vector<Thread *> m_vecThread;
};



#endif