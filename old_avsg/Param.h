#ifndef AVS_API_PARAM_H
#define AVS_API_PARAM_H

#include "Encap.h"
#include <unordered_set>
#include "BufPool.h"
#include "log/profile.h"
#include "log/progress_bar.h"

enum ArcType //压缩归档类型
{
    ARCTYPE_PE = 0,
    ARCTYPE_SE,
    ARCTYPE_MULTI
};

enum InFileType //文件类型
{
    INFILETYPE_FASTQ = 0,
    INFILETYPE_GZIP,
    INFILETYPE_MGZIP,
    INFILETYPE_PIPE
};

enum ActionType //执行动作
{
    ACTIONTYPE_DOINDEX = 0,
    ACTIONTYPE_DOENCODE,
    ACTIONTYPE_DODECODE,
    ACTIONTYPE_EXTRACT
};

enum AlignType //比对类型
{
    ALIGNTYPE_REFFREE = 0,
    ALIGNTYPE_BWA,
    ALIGNTYPE_HASH,
    ALIGNTYPE_MINIMIZER
};

enum OutFileType
{
    OUTFILETYPE_FILE = 0,
    OUTFILETYPE_PIPE
};

enum PARAM
{
    PARAM_DATACHECK = 1,
    PARAM_ALIGNTYPE,
    PARAM_ARCTYPE,
    PARAM_DEPENDENT,
    PARAM_IBLOCKNUM,
    PARAM_IMAXMIS,
    PARAM_IBLOCKSIZE,
    PARAM_KMER_SIZE,
    PARAM_SEQ_COMPRESS_MODE,
    PARAM_REFMD5,
    PARAM_FQFILE,
    PARAM_ADDPARAM,
    PARAM_BLOCKINFO
};

enum Extract
{
    EXTRACT_MATCH = 0,  //正则筛选
    EXTRACT_HEAD,       //提取前n条read
    EXTRACT_SRAND,      //随机提取n条read
    EXTRACT_LIST,       //显示压缩的多文件序号
    EXTRACT_GETFILES    //提取相应序号的文件
};

enum PipeOut
{
    PIPEOUT_NO = 0,
    PIPEOUT_SE,
    PIPEOUT_PE_1,
    PIPEOUT_PE_2,
    PIPEOUT_PE
};

enum Dependent
{
    DEPENDENT_NO = 0,    //无依赖
    DEPENDENT_SEQONQUAL, //seq依赖qual
    DEPENDENT_QUALONSEQ  //qual依赖seq
};

enum SeqCompressMode
{
    SEQ_PRESERVE_ORDER = 0,
    SEQ_REORDER
};

enum DataCheck
{
    DATACHECK_NO = 0, //不校验
    DATACHECK_MD5,
    DATACHECK_XXHASH
};

class Param
{
public:
    static Param *GetInstance()
    {
        static Param param;
        return &param;
    }
    void add(const string &str);
    int ParamEncap(char *buf);
    void ParamUnEncap(uint32_t size, char *buf);
private:
    Param();
    ~Param();
    void getRefMd5();
    void decompressParam(int id, int size, char *pbuf);
    int EncapRefMd5(uint32_t id, char *pbuf);
    int EncapBlockInfoSE(char *pbuf);
    int EncapBlockInfoPE(char *pbuf);
    int EncapBlockInfo(uint32_t id, char *pbuf);

    void addSEinfo(BlockInfoSE &info);
    int UnEncapBlockInfo(char *pbuf);
    int UnEncapBlockInfoPE(char *pbuf);
    int UnEncapBlockInfoSE(char *pbuf);
public:
    const uint32_t m_MinAlignLen = 64;
    const uint32_t m_MaxAlignLen = 500;
    // const uint32_t m_MaxAlignLen = 64;

    ActionType m_actionType;
    AlignType m_alignType = ALIGNTYPE_HASH;
    ArcType m_arcType;
    InFileType m_infiletype;
    OutFileType m_outfiletype;
    Dependent m_dependent = DEPENDENT_NO; //依赖
    DataCheck m_checktype = DATACHECK_MD5;
    uint32_t  m_iThreadnum = 1; //线程个数
    uint32_t  m_iBlocknum = 0; //压缩文件中包含的block个数
    uint32_t  m_iMaxmis = 8; //比对允许的最大mis个数
    uint32_t  m_iBlockSize = 250; //单位是M
    uint32_t  m_kmer_size  = 14; // kmer 大小,默认为14(在SRR554369中测试较好)
    SeqCompressMode m_seq_compress_mode = SEQ_PRESERVE_ORDER; // 碱基序列压缩模式
    PipeOut  m_PipeOut = PIPEOUT_NO; //管道输出
    std::string m_strRefFile; //参考序列路径
    std::vector<std::string> m_vecFqFile; //输入文件集合
    std::string m_strArcFile; //输出文件路径
    MemBufPool *m_memBufPoolPtr = nullptr; //内存池指针
    std::vector<BlockInfoSE> m_vecBlockInfoSE;
    std::vector<BlockInfoPE> m_vecBlockInfoPE;
    std::unordered_set<int> m_hashMultiple; //保存欲解压的多文件的序号
    std::vector<std::string> m_vecAddParam; //保存附加的参数
    std::vector<BlockMetaInfo> m_vecBlockMeta;
    uint64_t m_totalfilesize; 
    int m_outfd[2]; //输出文件句柄
    //std::string m_ref_filename; // 参考文件需要用户指定，不需要存储

    //筛选参数
    Extract m_extracttype;
    std::string m_strPattern; //正则表达式
    std::vector<int> m_vecReadCnt; //每个block需要提取的条数
private:
    unsigned char m_refmd5[16]; //保存参考序列的md5
    
};



#endif //AVS_API_PARAM_H