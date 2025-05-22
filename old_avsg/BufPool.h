#ifndef SEQARC_BUFPOOL_H
#define SEQARC_BUFPOOL_H

#include <vector>
#include <string>
#include <list>
#include <mutex>
#include <condition_variable>
#include "Memory.h"

class MemBuf
{
public:
    MemBuf(bool binit);
    ~MemBuf();
public:
    bool m_bonech = false; //第三行是否是单字符
    int m_penum = 0;       //pe读取的read个数
    int m_fileidx = 0;  //保存压缩文件的序号
    int m_blockidx = 0; //保存block块序号
    uint32_t m_len[2] = {0}; //记录block块大小
    uint64_t m_offset[2] = {0}; //记录block块偏移量
    char *m_bufptr[2] = {nullptr};
    Memory<int> m_vecpos1; //s1换行符的位置
    Memory<int> m_vecpos2; //s2换行符的位置
};

class MemBufPool
{
public:
    MemBufPool(int limit);
    ~MemBufPool();

    MemBuf *getEmptybuf(bool binit = true);
    MemBuf *getFullbuf();
    MemBuf *getIdxFullbuf(int idx);  //获取指定序号的内存
    void addEmptybuf(MemBuf *ptr);
    void addFullbuf(MemBuf *ptr);
    int getFullbufNum();
    int getTotalbufNum(){return m_totalbuf;}

private:
    int m_limit = 0;
    int m_totalbuf = 0;

    std::list<MemBuf *> m_fullList;   //保存内容为满的内存指针
    std::list<MemBuf *> m_emptyList;  //保存内容为空的内存指针
    std::mutex m_mtx_empty;
    std::condition_variable m_cv_empty;

    std::mutex m_mtx_full;
    std::condition_variable m_cv_full;
};

#endif //SEQARC_BUFPOOL_H