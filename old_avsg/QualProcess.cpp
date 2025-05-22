#include "QualProcess.h"
#include "IWorker.h"

QualProcess::QualProcess(IWorker *ptr)
{
    m_seqElPtr = ptr->getseqElPtr();
    m_qualElPtr = ptr->getqualElPtr();
    m_seqLenElPtr = ptr->getseqLenElPtr();
    m_profile = ptr->m_profile;
}

QualProcess::~QualProcess()
{
}

static constexpr size_t BlockCompressSizeLimit = 1000 * 1024 * 1024;

int QualProcess::simple_compress(char *outptr) {
    auto timer1 = m_profile->getTimer();
    timer1->start("QualMd5Time");
    int szmd5 = calcMd5(outptr);
    outptr += szmd5;
    timer1->stop();

    timer1->start("QualTime");
    int szqual = docompressQual(outptr);
    timer1->stop();
    return szmd5+szqual;
}

int QualProcess::one_block_qual_compress(uint8_t * psrc, uint8_t * psrc_seq, uint32_t * len_arr,
                                         uint64_t block_size,uint64_t read_num, char * outptr) 
{
    uint64_t destlen = block_size;
    double cal_row = destlen/double(read_num);
    int real_row = len_arr[0];
    int ret;
    if(cal_row==real_row)
        ret = aco.aco_compress((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);
    else
        ret = aco.aco_compress_unalign((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);
    return destlen;
}

int QualProcess::block_compress(char *outptr) {
    auto timer1 = m_profile->getTimer();
    timer1->start("QualMd5Time");
    int szmd5 = 0, szqual = 0;
    if (m_paramPtr->m_checktype == DATACHECK_MD5) {
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), m_szmd5);
        memcpy(outptr, m_szmd5, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &m_hashval);
        memcpy(outptr, &m_hashval, 16);
    }
    szmd5 = 16;
    outptr += szmd5;
    timer1->stop();

    timer1->start("QualTime");
    size_t t_num = 2; // 暂时固定为t_num块
    const uint32_t total_read_num = m_seqLenElPtr->m_LenArry.size();  //输出存起来

    const uint32_t unit_read_num = total_read_num / t_num + 1 ;
    uint32_t process_read_num = 0;

    std::vector<uint8_t*> psrc(t_num), psrc_seq(t_num);//存储碱基和质量分数

    std::vector<uint32_t*> len_arr(t_num);
    std::vector<uint64_t> block_size(t_num);
    std::vector<uint64_t> block_read_num(t_num);
    std::vector<std::vector<char>> buffer(t_num);

    auto * qual_ptr = reinterpret_cast<uint8_t*>(m_qualElPtr->m_rawBuf.getdata());
    auto * seq_ptr = reinterpret_cast<uint8_t*>(m_seqElPtr->m_rawBuf.getdata());  //获取数据

    uint32_t * len_ptr = m_seqLenElPtr->m_LenArry.getdata();

    int t_id = 0;
    while (t_id < t_num) 
    {
        uint64_t size = 0;
        psrc[t_id] = qual_ptr;
        psrc_seq[t_id] = seq_ptr;
        len_arr[t_id] = len_ptr;

        block_read_num[t_id] = (process_read_num + unit_read_num >= total_read_num) ? total_read_num - process_read_num : unit_read_num;
        for (auto i = 0; i < block_read_num[t_id]; ++i) 
        {
            size += len_ptr[i];
        }
        block_size[t_id] = size;
        qual_ptr += size;
        seq_ptr += size;
        len_ptr += block_read_num[t_id];
        process_read_num += block_read_num[t_id];
        t_id++;
    }
    for (size_t i = 0; i < t_num; ++i) {
        buffer[i].resize(block_size[i] / 2);
    }
#pragma omp parallel for
    for (size_t i = 0; i < t_num; ++i) 
    {
        int size = one_block_qual_compress(psrc[i],psrc_seq[i],len_arr[i],block_size[i],block_read_num[i],buffer[i].data());
        buffer[i].resize(size);
    }
    for (size_t i = 0; i < t_num; ++i) {
        uint32_t size = buffer[i].size();
        std::memcpy(outptr, &size, 4);
        outptr += 4;
        std::memcpy(outptr, buffer[i].data(), size);
        outptr += size;
        szqual += (4 + size);
    }
    timer1->stop();
    return szmd5 + szqual;
}

int QualProcess::compress(char *outptr)
{
    if (!m_paramPtr->m_alignType && m_paramPtr->m_iBlockSize >= BlockCompressSizeLimit) {
        return block_compress(outptr);
    } else {
        return simple_compress(outptr);
    }
}

int QualProcess::calcMd5(char *outptr)
{
    int idlen = Encap::setID(BLOCKQUAL_MD5, outptr);
    outptr += idlen;
    Encap::setSize(16, 1, outptr);
    outptr += 1;
    if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), m_szmd5);
        memcpy(outptr, m_szmd5, 16); //把数据变为MD5，存储在编码序列里面，后面会有调用
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &m_hashval);
        memcpy(outptr, &m_hashval, 16);
    }
    return 17+idlen;
}

int QualProcess::docompressQual(char *outptr)
{
    int idlen = Encap::setID(BLOCKQUAL_QUAL, outptr); //把对应的id封装进去，用来识别
    outptr += idlen;
    char *psize = outptr;
    outptr += SIZENUM4;

    uint64_t destlen = m_qualElPtr->m_rawBuf.size();
    uint64_t read_num = m_seqLenElPtr->m_LenArry.size();

    uint32_t *len_arr = m_seqLenElPtr->m_LenArry.getdata();
    uint8_t *psrc = (uint8_t*)m_qualElPtr->m_rawBuf.getdata(); //指向质量分数的指针，因为在编码，有数据存在

    // DEBUG->mmm
    uint8_t *psrc_seq = (uint8_t*)m_seqElPtr->m_rawBuf.getdata();

    double cal_row = destlen/double(read_num);
    int real_row = m_seqLenElPtr->m_LenArry[0];
    int ret;
    if(cal_row==real_row)
        ret = aco.aco_compress((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);
    else
        ret = aco.aco_compress_unalign((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);

    Encap::setSize(destlen, SIZENUM4, psize);
    return destlen+idlen+SIZENUM4;
}

int QualProcess::simple_decompress(int inlen, char *inptr) {
    char *data = inptr;
    EncapInfo info;
    uint32_t readlen = 0;
    while (readlen < inlen)
    {
        Encap::getValue(data, info);//inptr是编码数据
        int sid = info.val;
        data += info.len;
        readlen += info.len;
        Encap::getValue(data, info);
        data += info.len;
        readlen += info.len; //根据不同id解码出对应的数据

        dodecompressQual(sid, info.val, data);

        data += info.val;
        readlen += info.val;
    }
    return 0;
}

int QualProcess::block_decompress(int inlen, char *inptr) 
{
    if (m_paramPtr->m_checktype == DATACHECK_MD5) 
    {
        std::memcpy(m_szmd5, inptr, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        memcpy(&m_hashval, inptr, 16);
    }
    inptr += 16;//跳过前导信息
    inlen -= 16;
    size_t t_num = 0;
    std::vector<char*> data_ptr;
    std::vector<uint32_t> data_len;
    while (inlen > 0) 
    {
        uint32_t len = 0;
        std::memcpy(&len, inptr, 4);
        inptr += 4;
        inlen -= 4;
        data_ptr.push_back(inptr);//分区间进行
        data_len.push_back(len);
        inptr += len;
        inlen -= len;
        t_num++;
    }
    if (inlen != 0) {
        arc_stdout(level::ERROR,"Qual decompress error");
        std::exit(-1);
    }

    const uint32_t total_read_num = m_seqLenElPtr->m_LenArry.size();
    const uint32_t unit_read_num = total_read_num / t_num + 1;

    uint32_t process_read_num = 0;
    std::vector<uint8_t*> pdest(t_num), psrc_seq(t_num);
    std::vector<uint32_t*> len_arr(t_num);
    std::vector<uint64_t> block_size(t_num);
    std::vector<uint64_t> block_read_num(t_num);
    auto * qual_ptr = reinterpret_cast<uint8_t*>(m_qualElPtr->m_rawBuf.getdata());
    auto * seq_ptr = reinterpret_cast<uint8_t*>(m_seqElPtr->m_rawBuf.getdata());
    uint32_t * len_ptr = m_seqLenElPtr->m_LenArry.getdata();

    int t_id = 0;
    while (t_id < t_num) 
    {
        uint64_t size = 0;
        pdest[t_id] = qual_ptr;
        psrc_seq[t_id] = seq_ptr;
        len_arr[t_id] = len_ptr;
        block_read_num[t_id] = (process_read_num + unit_read_num >= total_read_num) ? total_read_num - process_read_num : unit_read_num;
        for (auto i = 0; i < block_read_num[t_id]; ++i) {
            size += len_ptr[i];
        }
        block_size[t_id] = size;
        qual_ptr += size;
        seq_ptr += size;
        len_ptr += block_read_num[t_id];
        process_read_num += block_read_num[t_id];
        t_id++;
    }

#pragma omp parallel for
    for (int i = 0; i < t_num; ++i) 
    {
        uint64_t destlen = block_size[i];
        double cal_row = block_size[i] /double(block_read_num[i]);
        int real_row = len_arr[i][0];   //  分区进行压缩，只取第一行，所以其实全部文件都用同一种，只不过分区进行

        if(cal_row==real_row)
            aco.aco_decompress(pdest[i],&destlen,(uint8_t*)data_ptr[i], psrc_seq[i],len_arr[i],data_len[i],block_read_num[i]);
        else
            aco.aco_decompress_unalign(pdest[i],&destlen,(uint8_t*)data_ptr[i], psrc_seq[i],len_arr[i],data_len[i],block_read_num[i],block_size[i]);
    }
    m_qualElPtr->m_rawBuf.updatesize(m_seqElPtr->m_rawBuf.size());
    checkData();
    return 0;
}

int QualProcess::decompress(int inlen, char *inptr)
{
    if (!m_paramPtr->m_alignType && m_paramPtr->m_iBlockSize >= BlockCompressSizeLimit) {
        return block_decompress(inlen, inptr);
    } else {
        return simple_decompress(inlen, inptr);
    }
}

int QualProcess::dodecompressQual(int id, int size, char *data)
{
    switch (id)
    {
        case BLOCKQUAL_MD5:
            memcpy(m_szmd5, data, 16); //把解码后的数据传入，比对MD5后的数据好很多，不用那么大
            break;
        case BLOCKQUAL_QUAL:
            decodeQual(size, data);
            checkData();
            break;
        default:
            break;
    }
    return 0;
}

int QualProcess::decodeQual(int inlen, char *inptr) //inlen是数据长度
{
    uint64_t destlen = m_qualElPtr->m_rawBuf.size(); //初始为0
    /*FILE *f = fopen("test.txt", "w");
    char *ptr = m_seqElPtr->m_rawBuf.getdata();
    for(int i=0;i<m_seqLenElPtr->m_LenArry.size();i++)
    {
        fwrite(ptr, 1, m_seqLenElPtr->m_LenArry[i], f);
        fwrite("\n", 1, 1, f);
        ptr += m_seqLenElPtr->m_LenArry[i];
    }
    fclose(f);

    FILE *ff = fopen("test2.txt", "w");
    fprintf(ff,"%lld",destlen);
    fclose(ff);*/

    uint64_t read_num = m_seqLenElPtr->m_LenArry.size();
    uint8_t *pdest = (uint8_t*)m_qualElPtr->m_rawBuf.getdata();//获取对应数据，没有则开辟对应的数据缓冲区，在函数里面进行赋值获取数据
    uint8_t *psrc_seq = (uint8_t*)m_seqElPtr->m_rawBuf.getdata();//已经解码了就有数据，未解码就获取空间用于输出
    uint32_t *len_arr = m_seqLenElPtr->m_LenArry.getdata();//每行长度
    uint64_t size_file = m_seqElPtr->m_rawBuf.size();

    double cal_row = size_file/double(read_num);
    int real_row = m_seqLenElPtr->m_LenArry[0];
    int ret;
    if(cal_row==real_row)
        ret = aco.aco_decompress(pdest,&destlen,(uint8_t*)inptr, psrc_seq,len_arr,inlen,read_num);
    else
        ret = aco.aco_decompress_unalign(pdest,&destlen,(uint8_t*)inptr, psrc_seq,len_arr,inlen,read_num,size_file);

    m_qualElPtr->m_rawBuf.updatesize(destlen);
    return destlen;
}

bool QualProcess::checkData()
{
    if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        unsigned char szmd5[16] = {0};
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), szmd5); //解码完成，调用数据MD5变化

        if(memcmp(m_szmd5, szmd5, 16) != 0)
        {
            arc_stdout(level::ERROR, "qual md5 check fail");
        }
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        XXH128_hash_t res;
        bool ret = calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &res);
        if(!ret || XXH128_isEqual(res, m_hashval)==0)
        {
            arc_stdout(level::ERROR, "name md5 check fail");
        }
    }
    return true;
}