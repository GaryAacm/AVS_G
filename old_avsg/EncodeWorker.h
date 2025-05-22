#ifndef AVS_API_ENCODEWORKER_H
#define AVS_API_ENCODEWORKER_H

#include "IWorker.h"
#include "BufPool.h"
#include <fstream>
#include "AlignInfo.h"

typedef struct _tagDegeInfo
{
    uint8_t max_dege_qual;
    uint32_t degecount;
}DegeInfo;

class EncodeWorker:public IWorker
{
public:
    EncodeWorker(int num, IRef* refptr);
    virtual ~EncodeWorker();
    virtual void doJob();
    virtual void doTask(MemBuf *ptr) = 0;
    void decomposeAlignInfo(AlignInfo &info);

    void cutHTDege(char *pseq, char *pqual, int len, AlignInfo &info);
    void getDegeInfo(char *pseq, int len);
    void getDegeInfo(char *pseq, int len, AlignInfo &info);
    void getDegeInfo(char *pseq, char *pqual, int len, AlignInfo &info);
    void getDegeInfo(char *pseq, char *pqual, int len);
    void DegeInfoProcess(char *pseq, char *pqual, int len);

    int DoCompressData();
    int DoCompressPart(IProcess *ptr, int id, char *outptr);
    
    uint64_t m_LastRefpos = 0;
protected:
    MemBufPool *m_memBufPoolPtr = nullptr;
    Memory<char> m_outbuf;
    DegeInfo m_degeinfo;

private:
    void reorder_process(const std::vector<uint32_t>& order);
};




#endif //AVS_API_ENCODEWORKER_H