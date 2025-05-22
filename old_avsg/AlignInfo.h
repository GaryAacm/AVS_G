#ifndef AVS_API_ALIGNINFO_H
#define AVS_API_ALIGNINFO_H

#include "util.h"
#include "memory.h"
#include "Param.h"

typedef struct _taghashseed
{
	uint32_t fqseed; //保存匹配上的seed
	uint32_t seqpos; //保存匹配上的seed对应的seqpos
	uint64_t seednum; //参考序列中seed对应的个数
}hashseed;

class AlignInfo
{
public:
    AlignInfo();
    ~AlignInfo();
    void init(int seqlen);
    void clear();
    void setprofile(std::shared_ptr<Module> &ptr){m_profile = ptr;}
public:
    bool m_bRev;    //是否是反向比对
    uint32_t m_misnum;   //变异碱基的个数（包括简并碱基）
    int m_seqlen;   //真实比对碱基的长度（剪切后的长度）
    uint32_t h_degelen;  //头部简并碱基剪切的长度
    uint32_t m_limit = 0; //每条序列mis的上限
    uint32_t m_seedextendcount; //比对的次数
    uint64_t m_alignpos; //比对上的参考序列的pos
    Memory<uint8_t> m_seqarry; //碱基序列二进制数组
    Memory<uint8_t> m_refarry; //参考序列二进制数组
    Memory<uint8_t> m_revseqarry; //存储反向互补碱基
    Memory<uint32_t> m_cigarL; //变异碱基的位置
    Memory<uint8_t> m_cigarV; //变异碱基的变化
    Memory<uint32_t> m_seqbase; //16个碱基保存成一个uint数据
    Memory<uint32_t> m_fqseed;  //存储碱基的seed值 
    Memory<uint32_t> m_revseqbase; //存储反向互补碱基的2进制数据
    Memory<uint32_t> m_revfqseed;  //存储反向互补碱基的seed值 

    hashseed m_hashseed[MAX_ALIGN_SEED];
    std::shared_ptr<Module> m_profile;
    int m_step = 0;
private:
    int m_maxlen = 512;
    Param *m_paramPtr = nullptr;
};




#endif