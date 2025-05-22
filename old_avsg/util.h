#ifndef AVS_API_UTIL_H
#define AVS_API_UTIL_H

#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <sys/types.h>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>
#include <string>
#include <inttypes.h>
#include "Memory.h"
#include "xxhash/xxh3.h"


using namespace std;

/********************************宏定义**********************/
#define DECODE_SHORT(a) ((a)[0] + ((a)[1]<<8))
#define DECODE_INT(a) ((a)[0] + ((a)[1]<<8) + ((a)[2]<<16) + ((a)[3]<<24))
#define RELEASEARRYPTR(p) {if(p) {delete[] p;p=nullptr;}}
#define RELEASEPTR(p) {if(p) {delete p; p=nullptr;}}
#define MODEL_SIZE(num) ((8+(num)*3))

/********************************常量**********************/
const int MAJOR_VERS = 1;
const int MINOR_VERS = 3;
const int STAGE_VERS = 0;

const int MAX_THREAD_NUM = 100;
const int PATH_LEN = 1024;
const int QMAX = 95;
const uint32_t SEED_MAX_DEPTH = 100000;
const uint32_t PART_SEED_MAX_DEPTH = 620;
const uint32_t MAX_ALIGN_SEED = 4;
const uint32_t MAX_ALIGN_POSITION = 3000;
const int ARCHEADLEN = 16;
const int SLEVEL = 3;
const int QBITS = 12;
const int QBITS_LESS = 10;

const int basesHash[4][4] = {3,1,0,2,
                             1,3,0,2,
                             1,0,3,2,
                             2,1,0,3};

const uint8_t hash_std_table[256] = {
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  0, 12,  1, 14,  4,  4,  2, 11,  4,  4,  8,  4,  5,  4,  4,
        4,  4,  6,  9,  3,  4, 13, 10,  4,  7,  4,  4,  4,  4,  4,  4,
        4,  0, 12,  1, 14,  4,  4,  2, 11,  4,  4,  8,  4,  5,  4,  4,
        4,  4,  6,  9,  3,  4, 13, 10,  4,  7,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,
        4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4
};
const char fqzdec[] = "ACGTNMRYKSWHBVD";

const uint64_t g_mask[33] = {0, 0x3, 0xf, 0x3f, 0xff, 0x3ff, 0xfff, 0x3fff, 0xffff, 0x3ffff, 0xfffff, 0x3fffff, 0xffffff,
                             0x3ffffff, 0xfffffff, 0x3fffffff, 0xffffffff, 0x3ffffffff, 0xfffffffff, 0x3fffffffff,
                             0xffffffffff, 0x3ffffffffff, 0xfffffffffff, 0x3fffffffffff, 0xffffffffffff,
                             0x3ffffffffffff, 0xfffffffffffff, 0x3fffffffffffff, 0xffffffffffffff,
                             0x3ffffffffffffff, 0xfffffffffffffff, 0x3fffffffffffffff, 0xffffffffffffffff};
const char mapCharry[4][4]={
    'G','C','T','N',
    'G','A','T','N',
    'C','A','T','N',
    'G','C','A','N'
};
const char mapBarry[4][4]={
    2,1,3,4,
    2,0,3,4,
    1,0,3,4,
    2,1,0,4
};
/********************************枚举**************************/

/********************************struct**********************/
typedef struct _tagBlockMetaInfo
{
    uint32_t fileidx;       //多文件压缩时的文件序号
    uint32_t compressSize;  //压缩后大小
    uint32_t originalSize1;  //压缩前大小
    uint32_t originalSize2;  //压缩前大小
    uint64_t compressOffset;  //压缩文件的偏移量
    uint64_t originalOffset1;  //原文件1的偏移量
    uint64_t originalOffset2;  //原文件2的偏移量
}BlockMetaInfo;


typedef struct _tagBlockInfoSE
{
    uint32_t fileidx;       //多文件压缩时的文件序号
    uint32_t compressSize;  //压缩后大小
    uint32_t originalSize;  //压缩前大小
    uint64_t compressoffset;//压缩文件偏移量
    uint64_t originaloffset;//压缩前文件偏移量
}BlockInfoSE;

typedef struct _tagBlockInfoPE
{
    uint32_t compressSize;  //压缩后大小
    uint32_t originalSize[2];  //压缩前大小
    uint64_t compressoffset;//压缩文件偏移量
    uint64_t originaloffset[2];//压缩前文件偏移量
}BlockInfoPE;

typedef struct _tagBlockParamInfo
{
    uint8_t onech; //read第三行是否是单字符
    uint32_t writeidx; //写入时的序号
    uint32_t blockidx; //block序号
    uint32_t readcount; //block的read条数
}BlockParamInfo;

/********************************Function**********************/
uint64_t getFileSize(const char *path);
bool isFileExist(const char *path);
bool isValidSrc(const char *path);
bool isValidRef(const char *path);
char rev(char ch);
int getbitnum(uint64_t data);
int int2bit(uint64_t data, int limit, Memory<uint8_t> &arry);
int int2bit(uint64_t data, int limit, uint8_t *buf);
void strSplit(const std::string &str, const std::string &delimiters, std::vector<string> &vectocken);
void splitPath(const char *path, string &dir, string &name);
void IntTo8Ch(uint64_t num, char *pbuf);
void IntTo4Ch(int num, char *pbuf);
void IntTo2Ch(int num, char *pbuf);
char revch(char ch);
int mapvar2base(int ch, int id);
void getFileName(string &str);
bool writeData(int fd, char *pbuf, uint64_t offset, uint64_t len);
bool calcDataHash(const void* buf, size_t len, XXH128_hash_t *result);
bool calcFileHash(const char *path, XXH128_hash_t *result);
#endif //AVS_API_UTIL_H