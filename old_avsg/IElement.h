#ifndef AVS_API_IELEMENT_H
#define AVS_API_IELEMENT_H

#include "Memory.h"
#include "Param.h"

class IElement
{
public:
    IElement()
    {
        m_paramPtr = Param::GetInstance();
    }
    virtual ~IElement(){}
    virtual void init() = 0;
    virtual void clear() = 0;
protected:

    Param *m_paramPtr = nullptr;
};



#endif //AVS_API_IELEMENT_H