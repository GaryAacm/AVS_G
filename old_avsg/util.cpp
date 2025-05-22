#include "util.h"
#include <stdarg.h>

uint64_t getFileSize(const char *path)
{
    uint64_t flength = 0;
    struct stat statbuf;
    int ret = stat(path, &statbuf);
    if (ret == -1) {
        return 0;
    }

    flength = statbuf.st_size;
    return flength;
}

bool isFileExist(const char *path)
{
    struct stat statbuf;
    int ret = stat(path, &statbuf);
    if (ret == -1) {
        return false;
    }
    return true;
}

bool isValidSrc(const char *path)
{
    string strRef(path);
    string strSuffix = strRef.substr(strRef.rfind('.')+1);
    if(strSuffix == "fq" || strSuffix == "fastq" || strSuffix == "gz")
    {
        return true;
    }
    return false;
}

bool isValidRef(const char *path)
{
    string strRef(path);
    string strSuffix = strRef.substr(strRef.rfind('.')+1);
    if(strSuffix == "fa" || strSuffix == "fasta")
    {
        return true;
    }
    return false;
}

char rev(char ch)
{
    switch(ch){
        case 'A':
            return 'T';
        case 'a':
            return 'T';
        case 'C':
            return 'G';
        case 'c':
            return 'G';
        case 'G':
            return 'C';
        case 'g':
            return 'C';
        case 'T':
            return 'A';
        case 't':
            return 'A';
        default:
            return ch;
    }
}

int getbitnum(uint64_t data)
{
    int count = 0;
    while (data > 0)
    {
        data >>= 1;
        count++;
    }
    return count;
}

int int2bit(uint64_t data, int limit, Memory<uint8_t> &arry)
{
    int count = 0;
    while (data > 0) {
        if (data & 1) {
            arry.push_back(1);
        } else {
            arry.push_back(0);
        }

        data >>= 1;
        count++;
    }

    if(limit < count)
    {
        printf("%s %d err data:%ld limit:%d\n", __FUNCTION__, __LINE__, data, limit);
        exit(1);
    }
    for (int i = count; i < limit; i++) {
        arry.push_back(0);
    }
    return 0;
}

int int2bit(uint64_t data, int limit, uint8_t *buf)
{
    int count = 0;
    while (data > 0) {
        if (data & 1) {
            buf[count] = 1;
        } else {
            buf[count] = 0;
        }

        data >>= 1;
        count++;
    }

    if(limit < count)
    {
        printf("%s %d err data:%ld limit:%d\n", __FUNCTION__, __LINE__, data, limit);
        exit(1);
    }
    for (int i = count; i < limit; i++) {
        buf[i] = 0;
    }
    return limit;
}

void strSplit(const std::string &str, const std::string &delimiters, std::vector<string> &vectocken)
{
    vectocken.clear();
    string::size_type lastpos = str.find_first_not_of(delimiters, 0);
    string::size_type pos = str.find_first_of(delimiters, lastpos);
    while (string::npos != lastpos || string::npos != pos)
    {
        vectocken.emplace_back(str.substr(lastpos, pos-lastpos));
        lastpos = str.find_first_not_of(delimiters, pos);
        pos = str.find_first_of(delimiters, lastpos);
    }
}

void splitPath(const char *path, string &dir, string &name)
{
    dir.clear();
    name.clear();
    string strpath(path);
    int pos = strpath.rfind('/');
    if(pos!=string::npos)
    {
        dir.append(strpath.substr(0, pos+1));
        name.append(strpath.substr(pos+1, strpath.length() -1));
    }
    else
    {
        name.append(path);
    }
}

void IntTo8Ch(uint64_t num, char *pbuf)
{
    *pbuf++ = (num >> 0) & 0xff;
    *pbuf++ = (num >> 8) & 0xff;
    *pbuf++ = (num >> 16) & 0xff;
    *pbuf++ = (num >> 24) & 0xff;

    *pbuf++ = (num >> 32) & 0xff;
    *pbuf++ = (num >> 40) & 0xff;
    *pbuf++ = (num >> 48) & 0xff;
    *pbuf++ = (num >> 56) & 0xff;
}

void IntTo4Ch(int num, char *pbuf)
{
    *pbuf++ = (num >> 0) & 0xff;
    *pbuf++ = (num >> 8) & 0xff;
    *pbuf++ = (num >> 16) & 0xff;
    *pbuf++ = (num >> 24) & 0xff;
}

void IntTo2Ch(int num, char *pbuf)
{
    *pbuf++ = (num >> 0) & 0xff;
    *pbuf++ = (num >> 8) & 0xff;
}

//A:0 C:1 G:2 T:3 N:4
char revch(char ch)
{
    switch (ch)
    {
    case 'a':
    case 'A':
        return 'T';
    case 'c':
    case 'C':
        return 'G';
    case 'g':
    case 'G':
        return 'C';
    case 't':
    case 'T':
        return 'A';
    default:
        return 'N';
    }
}

int mapvar2base(int ch, int id) {
    switch (ch) {
        case 0: {
            switch (id) {
                case 0:
                    return 2;
                case 1:
                    return 1;
                case 2:
                    return 3;
                case 3:
                    return 4;
            }
        }
        case 1: {
            switch (id) {
                case 0:
                    return 2;
                case 1:
                    return 0;
                case 2:
                    return 3;
                case 3:
                    return 4;
            }
        }
        case 2: {
            switch (id) {
                case 0:
                    return 1;
                case 1:
                    return 0;
                case 2:
                    return 3;
                case 3:
                    return 4;
            }
        }
        case 3: {
            switch (id) {
                case 0:
                    return 2;
                case 1:
                    return 1;
                case 2:
                    return 0;
                case 3:
                    return 4;
            }
        }
    }

    return 0;
}

void getFileName(string &str)
{
    int pos1 = str.rfind('/');
    int pos2 = str.find('.', pos1);
    int pos3 = str.find('.', pos2+1);
    pos1++;
    if(pos3 == string::npos)
    {
        str = str.substr(pos1);
    }
    else
    {
        str = str.substr(pos1, pos3-pos1);
    }
}

bool writeData(int fd, char *pbuf, uint64_t offset, uint64_t len)
{
    int ret = lseek(fd, offset, SEEK_SET);
    if(ret == -1)
    {
        printf("lseek data error\n");
        exit(1);
    }
    ret = write(fd, pbuf, len);
    if(ret == -1)
    {
        printf("write data error\n");
        exit(1);
    }
    return true;
}
bool calcDataHash(const void* buf, size_t len, XXH128_hash_t *result)
{
    bool bsuccess = true;
    XXH3_state_t* state = XXH3_createState();
    if(state)
    {
        if(XXH3_128bits_reset(state) != XXH_ERROR)
        {
            if(XXH3_128bits_update(state, buf, len) != XXH_ERROR)
            {
                *result = XXH3_128bits_digest(state);
            }
            else
            {
                bsuccess = false;
            }
        }
        else
        {
            bsuccess = false;
        }
        XXH3_freeState(state);
    }
    else
    {
        bsuccess = false;
    }

    return bsuccess;
}

// XXH128_hash_t ret={0,0};
// bool r = calcFileHash(path, &ret);
// char buf[64]={0};
// char *p = buf;
// unsigned char *rp = (unsigned char*)&ret;
// for (int i = 15; i >= 0; i--)
// {
//     sprintf(p, "%02x", rp[i]);
//     p+=2;
// }
// printf("xxhahs: %s\n", buf);
bool calcFileHash(const char *path, XXH128_hash_t *result)
{
    bool bsuccess = true;
    XXH3_state_t* state = XXH3_createState();
    if(state)
    {
        if(XXH3_128bits_reset(state) != XXH_ERROR)
        {
            unsigned char buffer[4096]; 
            FILE *file = fopen(path, "r");
            int len = 0;
            if(file)
            {
                while ((len = fread (buffer, 1, 4096, file)))  
                {
                    if(XXH3_128bits_update(state, buffer, len) == XXH_ERROR)
                    {
                        bsuccess = false;
                        break;
                    }
                }
                *result = XXH3_128bits_digest(state);
                fclose(file);
            }
            else
            {
                bsuccess = false;
            }
        }
        else
        {
            bsuccess = false;
        }
        XXH3_freeState(state);
    }
    else
    {
        bsuccess = false;
    }

    return bsuccess;
}