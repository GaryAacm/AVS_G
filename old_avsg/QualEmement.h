#ifndef AVS_API_QUALEMEMENT_H
#define AVS_API_QUALEMEMENT_H

#include "IElement.h"

class QualEmement:public IElement
{
public:
    QualEmement(){init();}
    ~QualEmement(){}
    virtual void init()
    {
        m_rawBuf.mem_malloc(m_paramPtr->m_iBlockSize/2);
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



#endif