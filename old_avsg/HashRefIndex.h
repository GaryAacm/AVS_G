
#ifndef AVS_API_HASHREFINDEX_H
#define AVS_API_HASHREFINDEX_H

#include "Param.h"

// const uint64_t g_mask[33] = {0, 0x3, 0xf, 0x3f, 0xff, 0x3ff, 0xfff, 0x3fff, 0xffff, 0x3ffff, 0xfffff, 0x3fffff, 0xffffff, 
// 														0x3ffffff, 0xfffffff, 0x3fffffff, 0xffffffff, 0x3ffffffff, 0xfffffffff, 0x3fffffffff, 
// 														0xffffffffff, 0x3ffffffffff, 0xfffffffffff, 0x3fffffffffff, 0xffffffffffff,
// 														0x3ffffffffffff, 0xfffffffffffff, 0x3fffffffffffff, 0xffffffffffffff,
// 														0x3ffffffffffffff, 0xfffffffffffffff, 0x3fffffffffffffff, 0xffffffffffffffff};
const int SEEDLEN = 14; //k-mer 长度
const int SEEDLIMIT = 16; //
const int INTERVAL = 2; //参考序列seed间隔

class IHashRefIndex{
public:
    IHashRefIndex(uint64_t filesize):
    m_iRefFileSize(filesize)
    {
        m_paramPtr = Param::GetInstance();
        m_ilimit = pow(2, SEEDLIMIT);
        m_ishort = SEEDLEN - 1;
        m_iSeedlen = SEEDLEN;
        m_iMaxnseeds = g_mask[SEEDLEN]+1;
    };
    virtual ~IHashRefIndex(){};
    virtual bool createRefIndexShm() = 0;
    virtual bool loadRefIndexShm() = 0;
    virtual bool loadRefIndexFile() = 0;

    virtual void setSeqint(uint64_t Seqint, uint8_t chval) = 0;
    virtual void setEndSeqint(uint64_t Chnum, uint64_t Seqint) = 0;
    virtual void setSeednum(uint32_t seedseq) = 0;
    virtual void setSeedind(uint64_t mask) = 0;
    virtual void setSeedpos(uint32_t seedseq, uint64_t Chnum) = 0;
    virtual bool writeIndexFile() = 0;
    virtual uint64_t getSeedlen() = 0;
    virtual uint64_t getChnum() = 0;
    virtual uint64_t getSeednum(uint32_t val) = 0;
    virtual uint32_t getSeqint(uint64_t val) = 0;
    virtual uint64_t getSeedpos(uint32_t val, uint64_t index) = 0;

    char *createShm(uint64_t lsize);
    int getHalf(){return m_ihalf;}
    virtual void initMemory() = 0;
protected:
    bool m_bNeedFree = true;
    Param *m_paramPtr = nullptr;
    int m_ihalf = 16; //32bit的一半
    int m_ishort = 0; //Seedlen-1
    uint32_t  m_ilimit = 0;
    uint32_t  m_iSeedlen = 0;
    uint64_t  m_iMaxnseeds = 0;
    uint64_t  m_iRefFileSize  = 0;
    uint32_t* m_pSeqint = nullptr;
};

class HashRefIndex32:public IHashRefIndex{
public:
    HashRefIndex32(uint64_t filesize);
    ~HashRefIndex32();

    virtual bool createRefIndexShm();
    virtual bool loadRefIndexShm();
    virtual bool loadRefIndexFile();

    virtual void setSeqint(uint64_t Seqint, uint8_t chval);
    virtual void setEndSeqint(uint64_t Chnum, uint64_t Seqint);
    virtual void setSeednum(uint32_t seedseq);
    virtual void setSeedind(uint64_t mask);
    virtual void setSeedpos(uint32_t seedseq, uint64_t Chnum);
    virtual bool writeIndexFile();
    virtual uint64_t getSeedlen(){return m_iSeedlen;};
    virtual uint64_t getChnum(){return m_iChnum;};
    virtual uint64_t getSeednum(uint32_t val){return m_pSeednum[val];};
    virtual uint32_t getSeqint(uint64_t val){return m_pSeqint[val];};
    virtual uint64_t getSeedpos(uint32_t val, uint64_t index){return m_pSeedpos[m_pSeedind[val]+index];};
    virtual void initMemory();
private:
    uint32_t  m_iChnum = 0; //参考序列碱基个数
    uint32_t  m_iSeqint = 0; //参考序列碱基二进制化字节数
    uint32_t  m_iTotnpos = 0;//总位置数
	uint32_t* m_pSeednum = nullptr;
	uint32_t* m_pSeedind = nullptr;
	uint32_t* m_pSeedpos = nullptr; 

    uint32_t* m_pSeedindshift = nullptr;
};

class HashRefIndex64:public IHashRefIndex{
public:
    HashRefIndex64(uint64_t filesize);
    ~HashRefIndex64();
    virtual bool createRefIndexShm();
    virtual bool loadRefIndexShm();
    virtual bool loadRefIndexFile();

    virtual void setSeqint(uint64_t Seqint, uint8_t chval);
    virtual void setEndSeqint(uint64_t Chnum, uint64_t Seqint);
    virtual void setSeednum(uint32_t seedseq);
    virtual void setSeedind(uint64_t mask);
    virtual void setSeedpos(uint32_t seedseq, uint64_t Chnum);
    virtual bool writeIndexFile();
    virtual uint64_t getSeedlen(){return m_iSeedlen;};
    virtual uint64_t getChnum(){return m_iChnum;};
    virtual uint64_t getSeednum(uint32_t val){return m_pSeednum[val];};
    virtual uint32_t getSeqint(uint64_t val){return m_pSeqint[val];};
    virtual uint64_t getSeedpos(uint32_t val, uint64_t index){return m_pSeedpos[m_pSeedind[val]+index];};
    virtual void initMemory();
private:
    uint64_t  m_iChnum = 0; //参考序列碱基个数
    uint64_t  m_iSeqint = 0; //参考序列碱基二进制化字节数
    uint64_t  m_iTotnpos = 0;//总位置数
	uint64_t* m_pSeednum = nullptr;
	uint64_t* m_pSeedind = nullptr;
	uint64_t* m_pSeedpos = nullptr; 

    uint64_t* m_pSeedindshift = nullptr;
};
#endif //AVS_API_HASHREFINDEX_H