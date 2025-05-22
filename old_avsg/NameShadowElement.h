#ifndef AVS_API_NAMESHADOWELEMENT_H
#define AVS_API_NAMESHADOWELEMENT_H

#include "IElement.h"

class NameShadowElement:public IElement
{
public:
    NameShadowElement(){init();}
    ~NameShadowElement(){}

    virtual void init()
    {
        m_LenArry.mem_malloc(1024*1024);
        m_rawBuf.mem_malloc(m_paramPtr->m_iBlockSize/5);
    }
    virtual void clear()
    {
        m_LenArry.clear();
        m_rawBuf.clear();
    }
 
    Memory<uint32_t> m_LenArry;
    Memory<char> m_rawBuf;

};



#endif //AVS_API_NAMESHADOWELEMENT_H