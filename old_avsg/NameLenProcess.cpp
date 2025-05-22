#include "NameLenProcess.h"
#include "IWorker.h"

NameLenProcess::NameLenProcess(IWorker *ptr)
{
    m_nameLenElPtr = ptr->getnamelenElPtr();
}

NameLenProcess::~NameLenProcess()
{
}

int NameLenProcess::compress(char *outptr)
{
    return 0;
}

int NameLenProcess::decompress(int inlen, char *inptr)
{
    return 0;
}
