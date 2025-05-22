#ifndef AVS_API_COMBINECOMPRESS_H
#define AVS_API_COMBINECOMPRESS_H

#include "IProcess.h"
#include "COMTXTCompress.h"

enum COMBINE
{
    BLOCKCOMBINEQUAL_MD5 = 1,
    BLOCKCOMBINESEQ_MD5,
    BLOCKCOMBINE
};
class IWorker;
class SeqLenElement;
class QualEmement;
class SeqElement;

class CombineCompress:public IProcess
{
public:
    CombineCompress(IWorker *ptr);
    ~CombineCompress();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);

private:
    int calcMd5_QUAL(char *outptr);
    int calcMd5_SEQ(char *outptr); 
    bool checkData();
    int docompressCombine(char *outptr);
    int dodecompressCombine(int id, int size, char *data);
    int decodeCombine(int inlen, char *inptr);

    int simple_compress(char *outptr);
    int simple_decompress(int inlen, char *inptr);

    int block_compress(char *outptr);
    int block_decompress(int inlen, char *inptr);

    int one_block_combine_compress(uint8_t * psrc, uint8_t * psrc_seq, uint32_t * len_arr, uint64_t block_size, uint64_t read_num, char * outptr);

private:
    unsigned char m_szmd5_qual[16];
    unsigned char m_szmd5_seq[16];
    XXH128_hash_t m_hashval_qual;
    XXH128_hash_t m_hashval_seq;
    QualEmement *m_qualElPtr = nullptr;
    SeqElement *m_seqElPtr = nullptr;
    SeqLenElement *m_seqLenElPtr = nullptr;
    std::shared_ptr<Module> m_profile;
    COMTXTCompress comtxt;
    uint64_t totsize_file = 0;
    uint64_t totread_num = 0;
};

#endif