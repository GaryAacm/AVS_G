#ifndef AVS_API_SEQELEMENT_H
#define AVS_API_SEQELEMENT_H

#include "IElement.h"

class SeqElement:public IElement
{
public:
    SeqElement(){init();}
    virtual ~SeqElement(){
    }
    virtual void init()
    {
        m_rawBuf.mem_malloc(m_paramPtr->m_iBlockSize/2);
        m_bufPtr = m_rawBuf.getdata();

        if (m_paramPtr->m_alignType) {
            m_vecRefIdx.mem_malloc(1024 * 1024);
            m_vecBlockPos.mem_malloc(10 * 1024 * 1024);
            m_vecCigaL.mem_malloc(1024 * 1024);
            m_vecMis.mem_malloc(1024 * 1024);
            m_vecRev.mem_malloc(1024 * 1024);
            m_vecCigaV.mem_malloc(1024 * 1024);
            m_vecPERel.mem_malloc(1024 * 1024);
            m_vecDegeTip.mem_malloc(1024 * 1024);
            m_vecDegeCh.mem_malloc(1024 * 1024);
            m_vecDegeMaxQual.mem_malloc(1024 * 1024);
            m_vecNDegeCnt.mem_malloc(1024 * 1024);
            m_vecNDegePos.mem_malloc(1024 * 1024);

            m_vecPosEqual.mem_malloc(1024 * 1024);
            m_vecBlockPos_wj.mem_malloc(1024 * 1024);
        }
    }
    virtual void clear()
    {
        m_InsrNum = 0;
        m_alignCount = 0;
        m_rawBuf.clear();
        m_bufPtr = m_rawBuf.getdata();

        if (m_paramPtr->m_alignType) {
            m_vecRefIdx.clear();
            m_vecBlockPos.clear();
            m_vecCigaL.clear();
            m_vecMis.clear();
            m_vecRev.clear();
            m_vecCigaV.clear();
            m_vecPERel.clear();
            m_vecDegeTip.clear();
            m_vecDegeCh.clear();
            m_vecDegeMaxQual.clear();
            m_vecNDegeCnt.clear();
            m_vecNDegePos.clear();

            m_vecPosEqual.clear();
            m_vecBlockPos_wj.clear();
        }
    }

    int m_InsrNum = 0;
    int m_alignCount = 0;
    char *m_bufPtr = nullptr;
    Memory<char> m_rawBuf;

    Memory<uint8_t> m_vecRefIdx; //参考序列block idx
    Memory<uint8_t> m_vecBlockPos; //比对pos二进制
    Memory<uint8_t> m_vecCigaL; //变异二进制
    Memory<uint8_t> m_vecMis;
    Memory<uint8_t> m_vecRev;
    Memory<uint8_t> m_vecCigaV;
    Memory<uint8_t> m_vecPERel; //保存pe关系

    Memory<uint8_t> m_vecDegeTip; //是否存在简并碱基
    Memory<uint8_t> m_vecDegeCh; //存放简并碱基字符
    Memory<uint8_t> m_vecDegeMaxQual; //存放简并最大质量值
    Memory<uint32_t> m_vecNDegeCnt; //非简并碱基个数
    Memory<uint32_t> m_vecNDegePos; //非简并碱基间隔
    
    Memory<uint8_t> m_vecPosEqual; //该pos是否和上一个相等
    Memory<uint32_t> m_vecBlockPos_wj; //比对pos
};




#endif