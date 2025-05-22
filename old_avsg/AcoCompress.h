#ifndef ACO_COMPRESS_H
#define ACO_COMPRESS_H

class AcoCompress
{
public:
    virtual int aco_compress(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num);
    virtual int aco_compress_unalign(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num);
    virtual int aco_decompress(uint8_t *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,int inlen,uint64_t read_num);
    virtual int aco_decompress_unalign(uint8_t *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,int inlen,uint64_t read_num,uint64_t size_file);

private:

};
// #include "IProcess.h"

// class IWorker;
// class SeqLenElement;
// class QualEmement;

// class QualProcess:public IProcess
// {
// public:
//     QualProcess(IWorker *ptr);
//     ~QualProcess();
//     virtual int compress(char *outptr);
//     virtual int decompress(int inlen, char *inptr);

// private:
//     QualEmement *m_qualElPtr = nullptr;
//     SeqLenElement *m_seqLenElPtr = nullptr;
//     std::shared_ptr<Module> m_profile;
// };

#endif