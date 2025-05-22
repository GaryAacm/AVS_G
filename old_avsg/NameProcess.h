/*
 * @Author: zhaozijian
 * @Date: 2021-05-07 14:42:06
 * @LastEditors: Please set LastEditors
 * @LastEditTime: 2021-07-23 11:20:07
 * @Description: file content
 */
#ifndef AVS_API_NAMEPROCESSW_H
#define AVS_API_NAMEPROCESSW_H

#include "IProcess.h"
#include "clr.h"
#include "simple_model.h"
#include <sstream>

enum BLOCKNAMEW
{
    BLOCKNAMEW_MD5 = 1,
    BLOCKNAMEW_NAME,
    BLOCKNAMEW_FCNT
};

enum BINTYPEW
{
    BINTYPEW_SEPARATOR = 0, //分隔符
    BINTYPEW_ALPHA, //字母
    BINTYPEW_DIGIT, //数字
    BINTYPEW_ALNUM, //数字+字母
    BINTYPEW_LENVAL, //seq长度值
    BINTYPEW_PETAIL, //pe尾部是1 和 2
    BINTYPEW_SETAIL, //se尾部是1 或 2
    BINTYPEW_LENGTH, //bin 是length
    BINTYPEW_ZERONUM, //bin是0开头的数字
};

enum DIGITSUB
{
    DIGITSUB_ZERO = 0, //差值为0，相同
    DIGITSUB_ONE, //差值为1，递增
    DIGITSUB_PLUS,
    DIGITSUB_MINUS
};

typedef struct
{
    uint32_t wtype:8;
    uint32_t wlen:24;
    uint32_t data;
    char abuf[64];
}win;

class IWorker;
class NameElement;
class NameLenElement;
class SeqLenElement;

class NameProcess:public IProcess
{
public:
    NameProcess(IWorker *ptr);
    ~NameProcess();
    virtual int compress(char *outptr);
    virtual int decompress(int inlen, char *inptr);

private:
    int calcMd5(char *outptr);
    int docompressName(char *outptr);
    int compressW(char *outptr);
    void splitName(char *name, int len);
    void initEnModel(char *name, int len);
    bool compressname(char *name, int len);
    int encodeFcnt(char *outptr);

    int dodecompress(int id, int size, char *data);
    int decodeFcnt(int size, char *data);
    int decompressW(char *data);
    int initDeModel(char *name);
    int decompressname(char *name, int idx);
    bool checkData();
private:
    int m_count = 0;
    XXH128_hash_t m_hashval;
    unsigned char m_szmd5[16];
    NameElement *m_nameElPtr = nullptr;
    NameLenElement *m_nameLenElPtr = nullptr;
    SeqLenElement *m_seqLenELPtr = nullptr;
    std::shared_ptr<Module> m_profile;

    Memory<win> m_winArry;
    SIMPLE_MODEL<128> m_model_f; //保存block的格式
    SIMPLE_MODEL<62> *m_model_m = nullptr; //保存字母
    SIMPLE_MODEL<16> *m_model_d = nullptr; //保存数字
    SIMPLE_MODEL<4> *m_model_t = nullptr; //保存指示
    Memory<uint64_t> m_vecFormat; //block里面每个foramt的个数 32:len|32:cnt

    RangeCoder m_rct;

    char *m_bufptr = nullptr; //指向输出内存的指针
    uint32_t m_flen = 0; //每个format的长度

    std::string m_strtmp;
};


#endif 