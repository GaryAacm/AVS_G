#ifndef AVS_API_SEQLENELEMENT_H
#define AVS_API_SEQLENELEMENT_H

#include "IElement.h"

class SeqLenElement:public IElement
{
public:
    SeqLenElement(){init();}
    ~SeqLenElement(){}
    virtual void init()
    {
        m_LenArry.mem_malloc(1024*1024);
    }
    virtual void clear()
    {
        m_LenArry.clear();
    }

    Memory<uint32_t> m_LenArry;
};


#endif //AVS_API_SEQLENELEMENT_H