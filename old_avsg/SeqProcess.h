#ifndef AVS_API_SEQPROCESS_H
#define AVS_API_SEQPROCESS_H

#include "IProcess.h"
#include "Encap.h"
#include "clr.h"
#include "seqmodel.h"
#include "simple_model.h"
#include "IRef.h"

enum BLOCKSEQ
{
    BLOCKSEQ_MD5 = 1,
    BLOCKSEQ_REFIDX,
    BLOCKSEQ_DEGETIP,
    BLOCKSEQ_DEGECH,
    BLOCKSEQ_DEGEMAXQUAL,
    BLOCKSEQ_NDEGECNT,
    BLOCKSEQ_DEGEPOS,
    BLOCKSEQ_PEINSRNUM,
    BLOCKSEQ_PERELATION,
    BLOCKSEQ_MAPPOS,
    BLOCKSEQ_MAPMIS,
    BLOCKSEQ_MAPREV,
    BLOCKSEQ_MAPCIGAL,
    BLOCKSEQ_MAPCIGAV,
    BLOCKSEQ_UNMAPSEQ
};

class IWorker;
class SeqElement;
class SeqLenElement;
class QualEmement;

class SeqProcess:public IProcess
{
public:
    SeqProcess(IWorker *ptr);
    ~SeqProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);
private:
    int calcMd5(char *outptr);
    int docompressSeq(char *outptr);
    void encode_seq(RangeCoder *rc, char *seq, int len);
    int compressPEInsrNum(char *outptr);
    int compressUnmapSeq(char *outptr);
    int compressSimpleModel(int range, int id, int count, uint8_t *parry, char *outptr);
    int compressRefidx(char *outptr);

    int compressDegeTip(char *outptr);
    int compressDegeCh(char *outptr);
    int compressDegeMaxQual(char *outptr);
    int compressNDegeCnt(char *outptr);
    int compressNDegePos(char *outptr);

    int compressPERelation(char *outptr);
    int compressAlignInfo_Pos(char *outptr);
    int compressAlignInfo_Mis(char *outptr);
    int compressAlignInfo_Rev(char *outptr);
    int compressAlignInfo_CigaL(char *outptr);
    int compressAlignInfo_CigaV(char *outptr);

    int decompressSimpleModel(int range, char* data, Memory<uint8_t> &arry);
    int dodecompressSeq(int id, uint32_t size, char* data);
    int decompressRefidx(char* data);

    int decompressDegeTip(char *data);
    int decompressDegeCh(char *data);
    int decompressDegeMaxQual(char *data);
    int decompressNDegeCnt(char *data);
    int decompressNDegePos(char *data);

    int decompressPEInsrNum(char* data);
    int decompressPERelation(char* data);
    int decompressAlignInfo_Pos(char* data);
    int decompressAlignInfo_Mis(char* data);
    int decompressAlignInfo_Rev(char* data);
    int decompressAlignInfo_CigaL(char* data);
    int decompressAlignInfo_CigaV(char* data);
    int decompressUnmapSeq(char* data);
    void decode_seq(bool bdege, RangeCoder *rc, char *seq, char *qual, int len);
    uint64_t getRealPos(bool isread2, uint64_t order);
    uint64_t getOffsePos(int count);
    int getCigal(int count);
    bool checkData();
    void kModelInit();
    void kModelEncode(RangeCoder *rc, uint64_t val);
    uint64_t kModelDecode(RangeCoder *rc);
    void encode_seq_formodel(char *seq, int len);
    void AlignInfoToSeq(bool bdege, char *pseq, char *pqual, int seqlen, int limit, uint64_t realpos);
    int getVal(uint8_t *&ptr, int count);
    void decodeAlignInfo(int seqlen, int cigar_num, uint64_t realpos);
    //seq independent
    int docompressSeq_independent(char *outptr);
    int compressUnmapSeq_independent(char *outptr);

    void c_i_divide(char *seq_p, uint32_t len, uint32_t &i, SIMPLE_MODEL<2>& seq_indicate, RangeCoder *rc);
    void d_i_divide(char *seq_p, uint32_t len, uint32_t &i, SIMPLE_MODEL<2>& seq_indicate, RangeCoder *rc);

    void encode_seq_normal(RangeCoder *rc, char *seq, int len);
    void encode_seq_dege(RangeCoder *rc, char *seq, int len, SIMPLE_MODEL<2> &seq_indicate);
    int decompressUnmapSeq_independent(char* data);
    void decode_seq_normal(RangeCoder *rc, char *seq, int len);
    void decode_seq_dege(RangeCoder *rc, char *seq, int len, SIMPLE_MODEL<2> &seq_indicate);
    void AlignInfoToSeq_independent(bool bdege, char *pseq, int seqlen, int limit, uint64_t realpos);
private:
    int m_seq_size;
    int NS_MASK;

    uint8_t *m_posPtr = nullptr;
    uint8_t *m_cigalPtr = nullptr;
    uint8_t *m_misptr = nullptr;
    std::vector<int> m_cigar_l;
    std::vector<int> m_cigar_v;
    uint32_t m_mapidx = 0;
    uint32_t m_cigavidx = 0;
    uint32_t m_degeidx = 0;
    uint32_t m_peRelidx = 0;
    uint32_t m_readfileidx = 0; 
    uint32_t m_degechidx = 0;
    uint32_t m_ndegeposidx = 0;

    uint64_t m_realpos = 0;
    unsigned char m_szmd5[16];
	XXH128_hash_t m_hashval;

    std::vector<SIMPLE_MODEL<2> *> m_vecKmodel;
    SIMPLE_MODEL<64> m_lenModel;
    uint8_t m_bitarry[64];

    SEQ_MODEL<uint8_t> *m_model_seq = nullptr;
    IRef *m_refPtr = nullptr;

    char *m_pmodelbuf = nullptr;

    SeqElement *m_seqElPtr = nullptr;
    SeqLenElement *m_seqLenElPtr = nullptr;
    QualEmement *m_qualElPtr = nullptr;
    Memory<uint8_t> m_refseq;  //保存查询比对上的参考序列
    std::shared_ptr<Module> m_profile;
};




#endif //AVS_API_SEQPROCESS_H