#ifndef AVS_API_QUALPROCESS_H
#define AVS_API_QUALPROCESS_H

#include "IProcess.h"
#include "AcoCompress.h"

enum BLOCKQUAL
{
    BLOCKQUAL_MD5 = 1,
    BLOCKQUAL_QUAL
};
class IWorker;
class SeqLenElement;
class QualEmement;
class SeqElement;

class QualProcess:public IProcess
{
public:
    QualProcess(IWorker *ptr);
    ~QualProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);

private:
    int calcMd5(char *outptr);
    bool checkData();
    int docompressQual(char *outptr);
    int dodecompressQual(int id, int size, char *data);
    int decodeQual(int inlen, char *inptr);

    int simple_compress(char *outptr);
    int simple_decompress(int inlen, char *inptr);

    int block_compress(char *outptr);
    int block_decompress(int inlen, char *inptr);

    int one_block_qual_compress(uint8_t * psrc, uint8_t * psrc_seq, uint32_t * len_arr, uint64_t block_size, uint64_t read_num, char * outptr);

private:
    unsigned char m_szmd5[16];
    XXH128_hash_t m_hashval;
    QualEmement *m_qualElPtr = nullptr;
    SeqElement *m_seqElPtr = nullptr;
    SeqLenElement *m_seqLenElPtr = nullptr;
    std::shared_ptr<Module> m_profile;
    AcoCompress aco;
};

#endif