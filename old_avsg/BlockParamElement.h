#ifndef AVS_API_BLOCKPARAMELEMENT_H
#define AVS_API_BLOCKPARAMELEMENT_H

#include "IElement.h"

class BlockParamElement:public IElement
{
public:
    BlockParamElement(){}
    ~BlockParamElement(){}
    virtual void init()
    {

    }
    virtual void clear()
    {
        m_onech = 0;
        m_longseq = 0;
        m_blockidx = 0;
        m_readcount = 0;
    }

    uint8_t m_onech; //read第三行是否是单字符
    uint8_t m_longseq; //是否是三代序列
    uint32_t m_blockidx; //block序号
    uint32_t m_readcount; //block的read条数
};


#endif //AVS_API_BLOCKPARAMELEMENT_H