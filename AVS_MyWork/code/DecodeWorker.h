#ifndef AVS_API_DECODEWORKER_H
#define AVS_API_DECODEWORKER_H

#include "IWorker.h"

typedef struct
{
    int size;
    char *pdata;
}DecodeData;


class DecodeWorker:public IWorker
{
public:
    DecodeWorker(int num, IRef* refptr);
    virtual ~DecodeWorker();
    virtual void doJob();
    virtual uint32_t recoverDataNoComment() = 0;
    virtual uint32_t recoverDataWithComment() = 0;
    virtual void outPutData(uint32_t len) = 0;

protected:
    int copyReadDataNoComment(int i, char *strout);
    int copyReadDataWithComment(int i, char *strout);
    void readData(uint32_t idx);
    void decodeData();
    // void doRealDecode(int id, uint32_t size, char* data);
protected:
    int m_count = 0;
    int m_readfileidx = 0;
    int m_blocksize[2] = {0};
    uint64_t m_blockoffset[2] = {0};
    std::fstream m_file;
    Memory<char> m_inbuf;
    DecodeData m_dparam, m_dseqlen, m_dnamelen, m_dname, m_dseq, m_dqual,m_combine;
private:
    void pe_reorder_process();
};


#endif
