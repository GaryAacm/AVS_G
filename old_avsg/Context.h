#ifndef AVS_API_CONTEXT_H
#define AVS_API_CONTEXT_H


class Param;
class IRef;
class FastqRead;
class ArcFile;

class Context
{
public:
    Context();
    ~Context();
    void initEncodeContext();
    void initDecodeContext();
    void initExtractContext();
    bool doBuildIndex();
private:
    void createRefPtr();
    
    bool doLoadIndex();
    void doReadAndEncode();
    void doDecode();
    void doExtract();
	void pipeOutDecode();
private:
    Param *m_paramPtr = nullptr;
    IRef *m_refPtr = nullptr;
};




#endif