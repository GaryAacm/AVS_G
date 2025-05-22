#ifndef AVS_API_ARCFILE_H
#define AVS_API_ARCFILE_H

#include "Param.h"
#include <map>

enum ARCPART
{
    ARCPART_HEADER = 1,
    ARCPART_BLOCK,
    ARCPART_TAIL,
};

enum ARCTAIL
{
    ARCTAIL_PARAM = 1,
};

class ArcFile
{
public:
    ArcFile();
    ~ArcFile();
    uint64_t writeFileInfo();
    int readBlockParam(std::map<uint32_t, BlockParamInfo> &mapBParam);
    int readTail();
private:
    int readInfoEncap();
    int writeParam();
    int doReadTail(int id, uint32_t size, char* data);
    
private:
    uint64_t m_totalBlockSize = 0;
    uint64_t m_filesize = 0;
    uint64_t m_tailoffset = 0;
    fstream m_file;
    Param *m_paramPtr = nullptr;
    Memory<char> m_tailbuf;
    char *m_bufPtr = nullptr;
    char m_readbuf[128]={0};
};



#endif //AVS_API_ARCFILE_H