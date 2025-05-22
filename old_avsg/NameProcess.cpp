/*
 * @Author: zhaozijian
 * @Date: 2021-05-07 14:42:06
 * @LastEditors: Please set LastEditors
 * @LastEditTime: 2021-07-23 11:22:46
 * @Description: file content
 */
#include "NameProcess.h"
#include "IWorker.h"

int8_t wch_to_val[128] = {
 -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
 -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
 -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,
  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  -1,  -1,  -1,  -1,  -1,  -1,
 -1, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24,
25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35,  -1,  -1,  -1,  -1,  -1,
 -1, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61,  -1,  -1,  -1,  -1,  -1
};

char wval_to_ch[64] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

static int hexcnt(int x) 
{
    int i=0;   
    while (x)
    {
        x>>=4;
        i++;
    }
    return i;
}
///////////////////////////////////////////////////


NameProcess::NameProcess(IWorker *ptr)
{
    m_nameElPtr = ptr->getnameElPtr();
    m_nameLenElPtr = ptr->getnamelenElPtr();
    m_seqLenELPtr = ptr->getseqLenElPtr();
    m_profile = ptr->m_profile;

    m_winArry.mem_calloc(128);
    m_vecFormat.mem_calloc(32);
    m_model_m = new SIMPLE_MODEL<62>[512];
    m_model_d = new SIMPLE_MODEL<16>[11];
    m_model_t = new SIMPLE_MODEL<4>[32];
}

NameProcess::~NameProcess()
{
    RELEASEARRYPTR(m_model_m);
    RELEASEARRYPTR(m_model_d);
    RELEASEARRYPTR(m_model_t);
}

int NameProcess::calcMd5(char *outptr)
{    
    int idlen = Encap::setID(BLOCKNAMEW_MD5, outptr);
    outptr += idlen;
    Encap::setSize(16, 1, outptr);
    outptr += 1;

    if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size(), &m_hashval);
        memcpy(outptr, &m_hashval, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        MDString(m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size(), m_szmd5);
        memcpy(outptr, m_szmd5, 16);
    }

    return 17+idlen;
}

int NameProcess::compress(char *outptr)
{
    auto timer1 = m_profile->getTimer();
    timer1->start("NameMd5Time");
    int szmd5 = calcMd5(outptr);
    outptr += szmd5;
    timer1->stop();

    timer1->start("NameTime");
    int szid = docompressName(outptr);
    outptr += szid;
    timer1->stop();

    int szcnt = encodeFcnt(outptr);
    return szmd5+szid+szcnt;
}

int NameProcess::docompressName(char *outptr)
{
    int idlen = Encap::setID(BLOCKNAMEW_NAME, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;

    m_count = m_seqLenELPtr->m_LenArry.size();
    int size = compressW(outptr);
    Encap::setSize(size, SIZENUM4, psize);
    return size+idlen+SIZENUM4;
}

int NameProcess::decompress(int inlen, char *inptr)
{
    char *data = inptr;
    EncapInfo info;
    uint32_t readlen = 0;
    while (readlen < inlen)
    {
        Encap::getValue(data, info);
        int sid = info.val;
        data += info.len;
        readlen += info.len;
        Encap::getValue(data, info);
        data += info.len;
        readlen += info.len;

        dodecompress(sid, info.val, data); 

        data += info.val;
        readlen += info.val;
    }

    return 0;
}

int NameProcess::dodecompress(int id, int size, char *data)
{
    switch (id)
    {
    case BLOCKNAMEW_MD5:
        if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
        {
            memcpy(&m_hashval, data, 16);
        }
        else if(m_paramPtr->m_checktype == DATACHECK_MD5)
        {
            memcpy(m_szmd5, data, 16);
        }
        break;
    case BLOCKNAMEW_NAME:
        m_bufptr = data;
        break;
    case BLOCKNAMEW_FCNT:
        decodeFcnt(size, data);

        m_count = m_seqLenELPtr->m_LenArry.size();
        decompressW(m_bufptr);
        checkData();
        break;
    default:
        break;
    }
    return 0;
}

int NameProcess::decodeFcnt(int size, char *data)
{
    m_vecFormat.clear();
    int cnt = size/8;
    uint64_t tmp = 0;
    for(int i=0;i<cnt;i++)
    {
        memcpy(&tmp, data, 8);
        data += 8;
        m_vecFormat.push_back(tmp);
    }
    return size;
}

bool NameProcess::checkData()
{
    if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        XXH128_hash_t res;
        bool ret = calcDataHash(m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size(), &res);
        if(!ret || XXH128_isEqual(res, m_hashval)==0)
        {
            arc_stdout(level::ERROR, "name md5 check fail");
        }
    }
    else if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        unsigned char szmd5[16] = {0};
        MDString(m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size(), szmd5);

        if(memcmp(m_szmd5, szmd5, 16) != 0)
        {
            arc_stdout(level::ERROR, "qual md5 check fail");
        }
    }
    return true;
}

int NameProcess::compressW(char *outptr)
{
    m_rct.output(outptr);
    m_rct.StartEncode();
    m_vecFormat.clear();
    m_winArry.clear();
    char *name_p = m_nameElPtr->m_rawBuf.getdata();
    char *ptmp;
    bool ret = false;
    uint64_t i = 0, oldi = 0;
    //FILE *f = fopen("raw.name", "w");
    RangeData rdata;
    while (i < m_count)
    {
        m_rct.midget(&rdata);
        ret = compressname(name_p, m_nameLenElPtr->m_LenArry[i]);
        if(!ret)
        {
            if(i)
            {
                uint64_t tmp = ((i-oldi)<<32) | m_flen;
                m_vecFormat.push_back(tmp);
                oldi = i;
            }
            m_rct.midset(&rdata);
            initEnModel(name_p, m_nameLenElPtr->m_LenArry[i]);
        }
        // fwrite(name_p, 1, m_nameLenElPtr->m_LenArry[i], f);
        // fwrite("\n", 1, 1, f);
        name_p += m_nameLenElPtr->m_LenArry[i];
        i++;
    }
    //fclose(f);

    uint64_t tmp = ((i-oldi)<<32) | m_flen;
    m_vecFormat.push_back(tmp);

    m_rct.FinishEncode();
    int sz = m_rct.size_out();
    return sz;
}

int NameProcess::encodeFcnt(char *outptr)
{
    int idlen = Encap::setID(BLOCKNAMEW_FCNT, outptr);
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;

    int size = m_vecFormat.size()*8;
    memcpy(outptr, m_vecFormat.getdata(), size);

    Encap::setSize(size, SIZENUM4, psize);
    return size + idlen + SIZENUM4;
}

void NameProcess::splitName(char *name, int len)
{
    //printf("###%.*s\n", len, name);
    m_winArry.clear();
    int8_t val = 0;
    uint32_t flen = 0, wlen = 0; //flen符号长度，wlen被符号分割的窗口长度
    uint8_t wtype = 0; //窗口类型
    int idx = 0, m=0, j = 0;
    bool blength = false;
    for(int i=0;i<len;i++)
    {
        val = wch_to_val[name[i]];
        if(val < 0) //符号
        {
            if(wlen)
            {
                m_winArry[idx].wlen = wlen;
                memcpy(m_winArry[idx].abuf, &name[m], wlen);
                if(strncasecmp("length", &name[m], wlen) == 0)
                {
                    blength = true;
                    j = idx;
                    m_winArry[idx].wtype = 7;
                }
                else
                {
                    m_winArry[idx].wtype = wtype;
                    if(wtype == 2)
                    {
                        if(m_winArry[idx].abuf[0]=='0')
                        {
                            m_winArry[idx].wtype = 8;
                        }
                        m_winArry[idx].data = atoi(&name[m]);
                    }
                }
                
                idx++;
                m = i;
            }
            flen++;
            wlen = 0;
            wtype = 0;
        }
        else
        {
            if(flen)
            {
                m_winArry[idx].wtype = wtype;
                m_winArry[idx].wlen = flen;
                memcpy(m_winArry[idx].abuf, &name[m], flen);
                idx++;
                m = i;
            }
            flen = 0;
            wlen++;
            wtype |= val>9 ? 1:2;
        }
    }

    m_winArry[idx].wtype = wtype;
    m_winArry[idx].wlen = wlen;
    if(wtype == 2)
    {
        if(blength && idx==j+2)
        {
            m_winArry[idx].wtype = 4;
        }

        if(wlen==1)
        {
            if(m_paramPtr->m_arcType == ARCTYPE_PE && 
                name[len-1]=='1' && name[2*len-1]=='2')
            {
                m_winArry[idx].wtype = 5;
            }
            if(m_paramPtr->m_arcType == ARCTYPE_SE &&
            (name[len-1]=='1' || name[len-1]=='2'))
            {
                m_winArry[idx].wtype = 6;
                m_winArry[idx].abuf[0] = name[len-1];
            }
        }
    }
    m_winArry.updatesize(idx+1);
}

void NameProcess::initEnModel(char *name, int len)
{
    m_flen = len;
    splitName(name, len);

    m_model_f.init();
    for(int i=0;i<len;i++)
    {
        m_model_f.encodeSymbol(&m_rct, name[i]);
    }

    for(int i=0;i<512;i++)
    {
        m_model_m[i].init();
    }
    for(int i=0;i<11;i++)
    {
        m_model_d[i].init();
    }

    for(int i=0;i<32;i++)
    {
        m_model_t[i].init();
    }
}

bool NameProcess::compressname(char *name, int len)
{
    if(m_winArry.size() == 0) return false;
    int k = 0, m = 0;
    uint8_t b = 0;
    char databuf[64]={0};
    for(int i=0;i<m_winArry.size();i++)
    {
        switch (m_winArry[i].wtype)
        {
        case BINTYPEW_SEPARATOR:
            if(memcmp(&name[k], m_winArry[i].abuf, m_winArry[i].wlen) !=0)
            {
                //printf("e show err\n");
                return false;
            }   
            k += m_winArry[i].wlen;
            break;
        case BINTYPEW_DIGIT:
        case BINTYPEW_ZERONUM:
        {
            int l=0;
            while (isdigit(name[k])&&k<len)
            {
                databuf[l++] = name[k++];
            }
            databuf[l]='\0';
            int val = atoi(databuf);
            int sub = val - m_winArry[i].data;
            m_winArry[i].data = val;

            if(sub == 0)
            {
                m_model_t[i].encodeSymbol(&m_rct, 0);
            }
            else if(sub == 1)
            {
                m_model_t[i].encodeSymbol(&m_rct, 1);
            }
            else
            {
                if(sub >0)
                {
                    m_model_t[i].encodeSymbol(&m_rct, 2);
                }
                else if (sub <0)
                {
                    m_model_t[i].encodeSymbol(&m_rct, 3);
                    sub = abs(sub);
                }

                int dl = hexcnt(sub);
                m_model_d[10].encodeSymbol(&m_rct, dl);
                int p = 0;
                while (sub)
                {
                    uint8_t h = sub&0xf;
                    sub >>= 4;
                    m_model_d[p++].encodeSymbol(&m_rct, h);
                }
            }
        }
            break;
        case BINTYPEW_ALPHA:
        case BINTYPEW_ALNUM:
            if(memcmp(&name[k], m_winArry[i].abuf, m_winArry[i].wlen) == 0) //相同
            {
                m_model_t[i].encodeSymbol(&m_rct, 1);
                k+= m_winArry[i].wlen;
            }
            else
            {
                m_model_t[i].encodeSymbol(&m_rct, 0);
                for(int j=0;j<m_winArry[i].wlen;j++,k++,m++)
                {
                    if(wch_to_val[name[k]]>=0)
                    {
                        m_model_m[m].encodeSymbol(&m_rct, wch_to_val[name[k]]);
                    }
                    else
                    {
                        //printf("e show err\n");
                        return false;
                    }
                }
            }
            break;
        case BINTYPEW_LENVAL: //长度值
            while (isdigit(name[k])&&k<len)
            {
                k++;
            }
            break;
        case BINTYPEW_PETAIL: //末尾的数字,指示se pe
        case BINTYPEW_SETAIL:
            k++;
            break;
        case BINTYPEW_LENGTH:
            k += m_winArry[i].wlen;
            break;
        default:
            break;
        }
    }

    if(k!=len)
    {
        arc_stdout(level::ERROR, "%s %d err", __FUNCTION__, __LINE__);
    }
    return true;
}

int NameProcess::decompressW(char *data) 
{
    uint32_t idlen = 0;
    m_rct.input(data);
    m_rct.StartDecode();

    char *name_p = m_nameElPtr->m_rawBuf.getdata();
    char *ptmp = name_p;
    for(int i=0;i<m_vecFormat.size();i++)
    {
        m_flen = m_vecFormat[i];
        uint32_t tcnt = m_vecFormat[i]>>32;

        initDeModel(name_p);
        name_p += m_flen;
        m_nameLenElPtr->m_LenArry.push_back(m_flen);
        

        for(int j=1;j<tcnt;j++)
        {
            idlen = decompressname(name_p, j);
            name_p += idlen;
            m_nameLenElPtr->m_LenArry.push_back(idlen);
        }
    }

    int size = name_p - ptmp;
    m_nameElPtr->m_rawBuf.updatesize(size);
    return size;
}

int NameProcess::initDeModel(char *name)
{
    m_model_f.init();
    for(int i=0;i<m_flen;i++)
    {
        name[i] = m_model_f.decodeSymbol(&m_rct);
    } 
    splitName(name, m_flen);

    for(int i=0;i<512;i++)
    {
        m_model_m[i].init();
    }
    for(int i=0;i<11;i++)
    {
        m_model_d[i].init();
    }

    for(int i=0;i<32;i++)
    {
        m_model_t[i].init();
    }
    return 0;
}

int NameProcess::decompressname(char *name, int idx)
{
    int k = 0, m = 0;
    uint8_t b = 0;
    char buf[64];
    for(int i=0;i<m_winArry.size();i++)
    {
        switch (m_winArry[i].wtype)
        {
        case BINTYPEW_SEPARATOR:
            memcpy(&name[k], m_winArry[i].abuf, m_winArry[i].wlen);
            k+=m_winArry[i].wlen;
            break;
        case BINTYPEW_DIGIT:
        {
            int t = m_model_t[i].decodeSymbol(&m_rct);
            int val = 0, dl = 0;
            switch (t)
            {
            case 0:
            {
                int len = sprintf(buf, "%d", m_winArry[i].data);
                memcpy(&name[k], buf, len);
                k += len;
            }
                break;
            case 1:
            {
                val = m_winArry[i].data + 1;
                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                memcpy(&name[k], buf, len);
                k += len; 
            }
                break;
            case 2:
            {
                dl = m_model_d[10].decodeSymbol(&m_rct);
                for(int p=0;p<dl;p++)
                {
                    uint8_t h = m_model_d[p].decodeSymbol(&m_rct);
                    val |= h<<(4*p);
                }
                val += m_winArry[i].data;

                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                memcpy(&name[k], buf, len);
                k += len; 
            }
                break;
            case 3:
            {
                dl = m_model_d[10].decodeSymbol(&m_rct);
                for(int p=0;p<dl;p++)
                {
                    uint8_t h = m_model_d[p].decodeSymbol(&m_rct);
                    val |= h<<(4*p);
                }
                val = m_winArry[i].data - val;

                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                memcpy(&name[k], buf, len);
                k += len; 
            }
                break;
            default:
                break;
            }
        }
            break; 
        case BINTYPEW_ALPHA:  
        case BINTYPEW_ALNUM:
            if(m_model_t[i].decodeSymbol(&m_rct))
            {
                memcpy(&name[k], m_winArry[i].abuf, m_winArry[i].wlen);
                k += m_winArry[i].wlen;
            }
            else
            {
                for(int j=0;j<m_winArry[i].wlen;j++,k++,m++)
                {
                    b = m_model_m[m].decodeSymbol(&m_rct);
                    name[k] = wval_to_ch[b];
                }
            }
            break;
        case BINTYPEW_LENVAL:
        {
            int len = sprintf(buf, "%d", m_seqLenELPtr->m_LenArry[idx]);
            memcpy(&name[k], buf, len);
            k += len;  
            break;
        }
        case BINTYPEW_PETAIL:
            name[k++]= (idx&1) ? '2': '1';
            break;
        case BINTYPEW_SETAIL:
            name[k++]= m_winArry[i].abuf[0];
            break;
        case BINTYPEW_LENGTH:
            memcpy(&name[k], m_winArry[i].abuf, m_winArry[i].wlen);
            k+=m_winArry[i].wlen;
            break;
        case BINTYPEW_ZERONUM:
        {
            int t = m_model_t[i].decodeSymbol(&m_rct);
            int val = 0, dl = 0;
            switch (t)
            {
            case 0:
            {
                memcpy(&name[k], m_winArry[i].abuf, m_winArry[i].wlen);
                k += m_winArry[i].wlen;
            }
                break;
            case 1:
            {
                val = m_winArry[i].data + 1;
                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                m_strtmp.clear();
                m_strtmp.append('0', m_winArry[i].wlen-len);
                m_strtmp.append(buf, len);
                memcpy(&name[k], m_strtmp.c_str(), m_winArry[i].wlen);
                k += m_winArry[i].wlen; 
            }
                break;
            case 2:
            {
                dl = m_model_d[10].decodeSymbol(&m_rct);
                for(int p=0;p<dl;p++)
                {
                    uint8_t h = m_model_d[p].decodeSymbol(&m_rct);
                    val |= h<<(4*p);
                }
                val += m_winArry[i].data;

                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                m_strtmp.clear();
                m_strtmp.append('0', m_winArry[i].wlen-len);
                m_strtmp.append(buf, len);
                memcpy(&name[k], m_strtmp.c_str(), m_winArry[i].wlen);
                k += m_winArry[i].wlen; 
            }
                break;
            case 3:
            {
                dl = m_model_d[10].decodeSymbol(&m_rct);
                for(int p=0;p<dl;p++)
                {
                    uint8_t h = m_model_d[p].decodeSymbol(&m_rct);
                    val |= h<<(4*p);
                }
                val = m_winArry[i].data - val;

                m_winArry[i].data = val;
                int len = sprintf(buf, "%d", m_winArry[i].data);
                m_strtmp.clear();
                m_strtmp.append('0', m_winArry[i].wlen-len);
                m_strtmp.append(buf,len);
                memcpy(&name[k], m_strtmp.c_str(), m_winArry[i].wlen);
                k += m_winArry[i].wlen; 
            }
                break;
            default:
                break;
            }
        }
        break;
        default:
            break;
        }
    }
    return k;
}

