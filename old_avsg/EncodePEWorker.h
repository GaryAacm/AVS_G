#ifndef AVS_API_ENCODEPEWORKER_H
#define AVS_API_ENCODEPEWORKER_H

#include "EncodeWorker.h"

typedef struct _tagInsrInfo
{
    uint32_t n;
    uint32_t med;
    uint32_t left;
    uint32_t right;
}InsrInfo;

typedef struct _tagPEAlign
{
    uint8_t rev[2];
    uint32_t mis[2];
    uint32_t seqlen[2];
    uint64_t alignpos[2];
}PEAlign;

class EncodePEWorker:public EncodeWorker
{
public:
    EncodePEWorker(int num, IRef* refptr);
    ~EncodePEWorker();
    virtual void doTask(MemBuf *ptr);
    void getBlockRead(MemBuf *ptr);
    void doAlign();
    void addBlockPEInfo();
private:
    bool CaclInsertSize(std::vector<int> &vec_insr);
    int AlignInfoProcessPE(PEAlign &peinfo);
private:
    uint32_t m_insertSize = 0; //pos1与pos2的差值范围
    int m_InsrNum = 0; //表示m_insertSize需要的bit个数
    InsrInfo m_insrInfo;
    BlockInfoPE m_info;
    AlignInfo m_alignInfo1;
    AlignInfo m_alignInfo2;
    PEAlign m_pealign;
    std::vector<int> m_vecinsr;
    std::vector<PEAlign> m_vecpealign;
};



#endif
