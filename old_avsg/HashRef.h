#ifndef AVS_API_HASHREF_H
#define AVS_API_HASHREF_H

#include "IRef.h"
#include "HashRefIndex.h"
#include "AlignInfo.h"


// typedef struct aligner_args{
// 	int kmerlen;
// 	int max_mismatch;
// 	uint64_t max_reference_len;
// 	int exp_mismatch;
// 	bool disable_slice;
// 	bool misopt;
// } aliargs;


class HashRef:public IRef
{
public:
    HashRef();
    ~HashRef();
    virtual bool buildRefIndex();
    virtual bool loadRefIndex();
    virtual bool freeRefIndex();
    bool isRefIndexExist();

    virtual uint32_t doSEAlign(AlignInfo &info);
    virtual uint32_t doPEAlign(AlignInfo &info1, AlignInfo &info2);
    virtual int query(uint64_t alignPos, int seqlen, uint8_t *refseq);

private:
    uint32_t getHashAlignInfo(AlignInfo &info);
    int getHashSeeds(AlignInfo &info);
    bool findHashSeeds(Memory<uint32_t> &fqseed, int start, int end, int maxdepth, hashseed &hseed, bool bfirst);
    bool gaplessHashAlignPositions(uint64_t refPos, AlignInfo &info);
    bool gaplessHashAlign(int idx, AlignInfo &info);
    bool hashAligner(AlignInfo &info);
    bool hashAlignerShortPart(AlignInfo &info);
private:
    bool m_bslice = true; 
    int m_kmerlen = 0;
    uint32_t m_reflen = 0;
    IHashRefIndex *m_hashIndexPtr = nullptr;
};


#endif //AVS_API_HASHREF_H