#ifndef COMTXTCOMPRESS_H
#define COMTXTCOMPRESS_H
#include "IProcess.h"
#include "Encap.h"
#include "clr.h"
#include "seqmodel.h"
#include "simple_model.h"
#include "IRef.h"

class COMTXTCompress
{
public:
    virtual int COMTXT_compress(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num);
    virtual int COMTXT_compress_unalign(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num);

    virtual int COMTXT_decompress(uint8_t *dest_qual,uint8_t *dest_seq,uint64_t *destLen_qual,uint64_t *destLen_seq,const uint8_t *source, const uint32_t *len_arr,int inlen,uint64_t read_num);
    virtual int COMTXT_decompress_unalign(uint8_t *dest_qual,uint8_t *dest_seq,uint64_t *destLen_qual,uint64_t *destLen_seq,const uint8_t *source, const uint32_t *len_arr,int inlen,uint64_t read_num,uint64_t size_file);

};

#endif