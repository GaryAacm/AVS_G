#ifndef AVS_API_FASTQREAD_H
#define AVS_API_FASTQREAD_H

#include <vector>
#include <string>
#include "Memory.h"
#include "Param.h"

class MemBufPool;
class MemBuf;
class IRef;
class ThreadPool;
class ProgressBar;

class FastqRead
{
public:
    FastqRead(IRef *refPtr=nullptr);
    ~FastqRead();
    void ReadFastq();
    void ReadPEFastq();
    void ReadPipe();
    bool isFinish();
    int getProcessCount();
    void setProcessBar(ProgressBar *barPtr){m_barPtr = barPtr;}
    int getBlockNum(){return m_blocknum;}
private:
    int getEndPos(int start, char *pbuf);
    void cultbuf(MemBuf *ptr);
    void cultPEbuf(MemBuf *ptr, bool bread);
    void startEncodeThread();
    void sendEndMsg();
    void getFirstLine(char *pbuf);

private:
    int m_blocknum = 0; //读取的block总数
    bool m_bonech = true;
    std::string m_strfirstline;
    std::string m_leftstr[2];
    Memory<int> m_vecpos1;
    Memory<int> m_vecpos2;
    MemBufPool *m_memBufPoolPtr = nullptr;
    IRef *m_refPtr = nullptr;
    Param *m_paramPtr = nullptr;
    ThreadPool *m_thPool = nullptr;
    ProgressBar *m_barPtr = nullptr;
};




#endif //AVS_API_FASTQREAD_H