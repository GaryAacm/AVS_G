#ifndef AVS_API_NAMELENELEMENT_H
#define AVS_API_NAMELENELEMENT_H

#include "IElement.h"

class NameLenElement:public IElement
{
public:
    NameLenElement(){init();}
    virtual ~NameLenElement(){}

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



#endif //AVS_API_NAMELENELEMENT_H