#ifndef AVS_API_NAMEELEMENT_H
#define AVS_API_NAMEELEMENT_H

#include "IElement.h"


class NameElement:public IElement
{
public:
    NameElement(){init();}
    virtual ~NameElement(){};

    virtual void init()
    {
        m_rawBuf.mem_malloc(m_paramPtr->m_iBlockSize);
        m_bufPtr = m_rawBuf.getdata();
    }
    virtual void clear()
    {
        m_rawBuf.clear();
        m_bufPtr = m_rawBuf.getdata();
    }

    Memory<char> m_rawBuf; 
    char *m_bufPtr = nullptr;
};



#endif //AVS_API_NAMEELEMENT_H