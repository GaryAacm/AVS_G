#include "CombineCompress.h"
#include "IWorker.h"

CombineCompress::CombineCompress(IWorker *ptr) //初始化得到数据
{
    m_seqElPtr = ptr->getseqElPtr();
    m_qualElPtr = ptr->getqualElPtr();
    m_seqLenElPtr = ptr->getseqLenElPtr();
    m_profile = ptr->m_profile;
}

CombineCompress::~CombineCompress()
{
}

static constexpr size_t CombineBlockCompressSizeLimit = 1000 * 1024 * 1024;

int CombineCompress::compress(char *outptr)
{
    if (!m_paramPtr->m_alignType && m_paramPtr->m_iBlockSize >= CombineBlockCompressSizeLimit) {
        return block_compress(outptr);
    } else {
        return simple_compress(outptr);
    }
}

int CombineCompress::decompress(int inlen, char *inptr)
{
    if (!m_paramPtr->m_alignType && m_paramPtr->m_iBlockSize >= CombineBlockCompressSizeLimit) 
    {
        return block_decompress(inlen, inptr);
    } 
    else 
    {
        return simple_decompress(inlen, inptr);
    }
}

int CombineCompress::simple_compress(char *outptr) {
    auto timer1 = m_profile->getTimer();
    timer1->start("CombineMd5Time");
    int szmd51 = calcMd5_QUAL(outptr);
    outptr += szmd51;

    int szmd52 = calcMd5_SEQ(outptr);
    outptr += szmd52;

    timer1->stop();

    timer1->start("CombineTime");
    int szcombine = docompressCombine(outptr);
    timer1->stop();
    return szmd51+szcombine+szmd52;
}

//================== Md5 ======================//

int CombineCompress::calcMd5_QUAL(char *outptr)
{
    int idlen = Encap::setID(BLOCKCOMBINEQUAL_MD5, outptr);
    outptr += idlen;
    Encap::setSize(16, 1, outptr);
    outptr += 1;
    if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), m_szmd5_qual);
        memcpy(outptr, m_szmd5_qual, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &m_hashval_qual);
        memcpy(outptr, &m_hashval_qual, 16);
    }

    return 17+idlen;
}

int CombineCompress::calcMd5_SEQ(char *outptr)
{
    int idlen = Encap::setID(BLOCKCOMBINESEQ_MD5, outptr);
    outptr += idlen;
    Encap::setSize(16, 1, outptr);
    outptr += 1;
    if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), m_szmd5_seq);
        memcpy(outptr, m_szmd5_seq, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &m_hashval_seq);
        memcpy(outptr, &m_hashval_seq, 16);
    }

    return 17+idlen;
}




int CombineCompress::docompressCombine(char *outptr)
{
    int idlen = Encap::setID(BLOCKCOMBINE, outptr);  //id用来识别
    outptr += idlen;
    char *psize = outptr;  //记录压缩前的位置
    outptr += SIZENUM4;

    
//============================data Initialization===========================//

    uint64_t destlen = m_qualElPtr->m_rawBuf.size();
    uint64_t destlen_seq = m_seqElPtr->m_rawBuf.size();
    uint64_t read_num = m_seqLenElPtr->m_LenArry.size();
    uint32_t *len_arr = m_seqLenElPtr->m_LenArry.getdata();

    uint8_t *psrc_qual = (uint8_t*)m_qualElPtr->m_rawBuf.getdata();
    uint8_t *psrc_seq = (uint8_t*)m_seqElPtr->m_rawBuf.getdata();

    uint64_t size_file = m_seqElPtr->m_rawBuf.size();
    //totread_num = read_num ;
    //totsize_file = size_file ;

    FILE *f = fopen("Simple_compress_file","w");
    fprintf(f,"%ld",read_num);
    fwrite("\n", 1, 1, f);
    fprintf(f,"%ld",size_file);
    fwrite("\n", 1, 1, f);
    fclose(f);

    FILE *ff = fopen("len_arr","w");
    for(int i = 0;i < read_num;i++)
    {
        fprintf(ff,"%u",len_arr[i]);
        fwrite("\n", 1, 1, ff);
    }
    fclose(ff);
    


//============================== process ===============================================//

    double cal_row = destlen/double(read_num);
    int real_row = m_seqLenElPtr->m_LenArry[0];
    int ret;
    if(cal_row==real_row)
        ret = comtxt.COMTXT_compress((char*)outptr,&destlen,psrc_qual, psrc_seq,len_arr,read_num);
    else
        ret = comtxt.COMTXT_compress_unalign((char*)outptr,&destlen,psrc_qual, psrc_seq,len_arr,read_num);

    Encap::setSize(destlen, SIZENUM4, psize); //重新调整大小，根据destlen的变化
    return destlen+idlen+SIZENUM4;
}

bool CombineCompress::checkData() //在解码之后比对数据
{
    if(m_paramPtr->m_checktype == DATACHECK_MD5)
    {
        unsigned char szmd5_qual[16] = {0};
        unsigned char szmd5_seq[16] = {0};
       
        MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), szmd5_seq);
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), szmd5_qual);

        if(memcmp(m_szmd5_seq, szmd5_seq, 16) != 0)
        {
            arc_stdout(level::ERROR, "seq md5 check fail");
        }

        if(memcmp(m_szmd5_qual, szmd5_qual, 16) != 0)
        {
            arc_stdout(level::ERROR, "qual md5 check fail");
        }

    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        XXH128_hash_t res_qual;
         XXH128_hash_t res_seq;
        bool ret_seq = calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &res_seq);
        bool ret_qual = calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &res_qual);

        if(!ret_seq || XXH128_isEqual(res_seq, m_hashval_seq)==0)
        {
            arc_stdout(level::ERROR, "seq name md5 check fail");
        }

        if(!ret_qual || XXH128_isEqual(res_qual, m_hashval_qual)==0)
        {
            arc_stdout(level::ERROR, "qual name md5 check fail");
        }
    }

    return true;
}

//============================== 分区压缩 ============================== //

int CombineCompress::block_compress(char *outptr) 
{
    auto timer1 = m_profile->getTimer();
    timer1->start("CombineMd5Time");
    int szmd5 = 0, szcombine = 0;
    if (m_paramPtr->m_checktype == DATACHECK_MD5) {
        MDString(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), m_szmd5_qual);
        memcpy(outptr, m_szmd5_qual, 16);
        outptr +=16;
        MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), m_szmd5_seq);
        memcpy(outptr, m_szmd5_seq, 16);
        outptr +=16;
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        calcDataHash(m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size(), &m_hashval_qual);
        memcpy(outptr, &m_hashval_qual, 16);
        outptr +=16;
        calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &m_hashval_seq);
        memcpy(outptr, &m_hashval_seq, 16);      
        outptr +=16;
    }

    szmd5 = 16 + 16;
    //outptr += szmd5;
    timer1->stop();

    timer1->start("CombineTime");
    size_t t_num = 2; // 暂时固定为t_num块
    const uint32_t total_read_num = m_seqLenElPtr->m_LenArry.size();
    const uint32_t unit_read_num = total_read_num / t_num + 1 ;
    uint64_t size_file = m_seqElPtr->m_rawBuf.size();
    uint32_t process_read_num = 0;
    std::vector<uint8_t*> psrc(t_num), psrc_seq(t_num);
    std::vector<uint32_t*> len_arr(t_num);
    std::vector<uint64_t> block_size(t_num);
    std::vector<uint64_t> block_read_num(t_num);
    std::vector<std::vector<char>> buffer(t_num);
    auto * qual_ptr = reinterpret_cast<uint8_t*>(m_qualElPtr->m_rawBuf.getdata());
    auto * seq_ptr = reinterpret_cast<uint8_t*>(m_seqElPtr->m_rawBuf.getdata());

    uint32_t * len_ptr = m_seqLenElPtr->m_LenArry.getdata();
    int t_id = 0;

    FILE *f = fopen("Block_len_arr.txt","w");
    for(int i=0;i<total_read_num;i++)
    {
        fprintf(f,"%u",len_ptr[i]);
        fwrite("\n", 1, 1, f);

    }
    fclose(f);

    FILE *ff = fopen("Block_compress_file.txt","w");
    fprintf(ff,"%u",total_read_num);
    fwrite("\n", 1, 1, ff);
    fprintf(ff,"%u",size_file);
    fwrite("\n", 1, 1, ff);
    fclose(ff);

    while (t_id < t_num) 
    {
        uint64_t size = 0;
        psrc[t_id] = qual_ptr;
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
    for (size_t i = 0; i < t_num; ++i) {
        buffer[i].resize(block_size[i] / 2);
    }
    //======================分好区间=====================//

#pragma omp parallel for
    for (size_t i = 0; i < t_num; ++i) 
    {
        int size = one_block_combine_compress(psrc[i],psrc_seq[i],len_arr[i],block_size[i],block_read_num[i],buffer[i].data());
        buffer[i].resize(size);
    }
    for (size_t i = 0; i < t_num; ++i) 
    {
        uint32_t size = buffer[i].size();
        std::memcpy(outptr, &size, 4);
        outptr += 4;
        std::memcpy(outptr, buffer[i].data(), size);
        outptr += size;
        szcombine += (4 + size);
    }

    timer1->stop();

    return szmd5 + szcombine;
}


int CombineCompress::one_block_combine_compress(uint8_t *psrc, uint8_t *psrc_seq, uint32_t *len_arr, uint64_t block_size, uint64_t read_num, char *outptr)
{
    uint64_t destlen = block_size;
    double cal_row = destlen/double(read_num);
    int real_row = len_arr[0];
    int ret;
    if(cal_row==real_row)
        ret = comtxt.COMTXT_compress((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);
    else
        ret = comtxt.COMTXT_compress_unalign((char*)outptr,&destlen,psrc, psrc_seq,len_arr,read_num);
    return destlen;

}


//======================== decompress ==============================//

int CombineCompress::simple_decompress(int inlen, char *inptr) 
{
    char *data = inptr;
    EncapInfo info;
    uint32_t readlen = 0;
    //FILE *fs = fopen("decompress_num.txt","w");
    //fprintf(fs,"%d 555\n",inlen);
    while (readlen < inlen)
    {
        Encap::getValue(data, info);//inptr是编码数据
        int sid = info.val;
        data += info.len;
        readlen += info.len;
        Encap::getValue(data, info);
        data += info.len;
        readlen += info.len; //根据不同id解码出对应的数据
        
        dodecompressCombine(sid, info.val, data);

        data += info.val;
        readlen += info.val;
    }
    //fclose(fs);
    return 0;
}


int CombineCompress::dodecompressCombine(int id, int size, char *data) //根据不同id解码对应的数据
{
   FILE *fss = fopen("decompress_id_use.txt","w");
   int num=8;
    fprintf(fss,"%d \n",num);
    switch (id)
    {
        case BLOCKCOMBINEQUAL_MD5:
            memcpy(m_szmd5_qual, data, 16); //把解码后的数据传入，比对MD5后的数据好很多，不用那么大
            break;
        case BLOCKCOMBINESEQ_MD5:
            memcpy(m_szmd5_seq, data, 16); //把解码后的数据传入，比对MD5后的数据好很多，不用那么大
            break;
        case BLOCKCOMBINE:
            decodeCombine(size, data);
            checkData();
            break;
        default:
            break;
    }
    fclose(fss);
    return 0;
}

int CombineCompress::decodeCombine(int inlen, char *inptr)  //inlen是数据长度 ,对下面数据进行文件读取
{
    uint64_t destlen_qual = m_qualElPtr->m_rawBuf.size(); //初始为0
    uint64_t destlen_seq = m_seqElPtr->m_rawBuf.size();

    uint8_t *dest_qual = (uint8_t*)m_qualElPtr->m_rawBuf.getdata();//获取对应数据，没有则开辟对应的数据缓冲区，在函数里面进行赋值获取数据
    uint8_t *dest_seq = (uint8_t*)m_seqElPtr->m_rawBuf.getdata();//已经解码了就有数据，未解码就获取空间用于输出

    uint64_t read_num  ; //= m_seqLenElPtr->m_LenArry.size();
    uint64_t size_file ; //= m_seqElPtr->m_rawBuf.size();

    FILE *f = fopen("Simple_compress_file","r");
    fscanf(f,"%ld\n%ld",&read_num,&size_file);
    fclose(f);

    uint32_t data[read_num];
    uint32_t *len_arr = data ;//每行长度
    FILE *ff = fopen("len_arr","r");
    for(int i=0;i<read_num;i++)
    {
        fscanf(ff,"%u",len_arr+i);
    }
    fclose(ff);
    
    


    double cal_row = size_file/double(read_num);
    int real_row = m_seqLenElPtr->m_LenArry[0];
    int ret;
    
    if(cal_row==real_row)
        ret = comtxt.COMTXT_decompress(dest_qual,dest_seq,&destlen_qual,&destlen_seq,(uint8_t*)inptr,len_arr,inlen,read_num);
    else
        ret = comtxt.COMTXT_decompress_unalign(dest_qual,dest_seq,&destlen_qual,&destlen_seq,(uint8_t*)inptr, len_arr,inlen,read_num,size_file);

    m_qualElPtr->m_rawBuf.updatesize(destlen_qual);
    m_seqElPtr->m_rawBuf.updatesize(destlen_seq);

    return destlen_qual + destlen_seq;
}

int CombineCompress::block_decompress(int inlen, char *inptr)
{
    if (m_paramPtr->m_checktype == DATACHECK_MD5) 
    {
        std::memcpy(m_szmd5_qual, inptr, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        memcpy(&m_hashval_qual, inptr, 16);
    }
    inptr += 16;//跳过前导信息
    inlen -= 16;

    if (m_paramPtr->m_checktype == DATACHECK_MD5) 
    {
        std::memcpy(m_szmd5_seq, inptr, 16);
    }
    else if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
    {
        memcpy(&m_hashval_seq, inptr, 16);
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
        arc_stdout(level::ERROR,"Combine decompress error");
        std::exit(-1);
    }

    uint32_t process_read_num = 0;
    std::vector<uint8_t*> pdest_qual(t_num), psrc_seq(t_num);
    std::vector<uint32_t*> len_arr(t_num);
    std::vector<uint64_t> block_size(t_num);
    std::vector<uint64_t> block_read_num(t_num);

    auto * qual_ptr = reinterpret_cast<uint8_t*>(m_qualElPtr->m_rawBuf.getdata());
    auto * seq_ptr = reinterpret_cast<uint8_t*>(m_seqElPtr->m_rawBuf.getdata());

    //uint32_t * len_ptr ; //= m_seqLenElPtr->m_LenArry.getdata();
    uint32_t total_read_num ; // = m_seqLenElPtr->m_LenArry.size();

    FILE *f = fopen("Block_compress_file","r");
    fscanf(f,"%ld\n",&total_read_num);
    fclose(f);

    uint32_t data[total_read_num];
    uint32_t *len_ptr = data ;//每行长度
    FILE *ff = fopen("Block_len_arr","r");
    for(int i=0;i<total_read_num;i++)
    {
        fscanf(ff,"%u",len_ptr+i);
    }
    fclose(ff);


    int t_id = 0;
    const uint32_t unit_read_num = total_read_num / t_num + 1;
    while (t_id < t_num) 
    {
        uint64_t size = 0;
        pdest_qual[t_id] = qual_ptr;
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
        uint64_t destlen_qual = block_size[i];
        uint64_t destlen_seq = block_size[i];
        double cal_row = block_size[i] /double(block_read_num[i]);
        int real_row = len_arr[i][0];   //  分区进行压缩，只取第一行，所以其实全部文件都用同一种，只不过分区进行

        if(cal_row==real_row)
            comtxt.COMTXT_decompress(pdest_qual[i],psrc_seq[i],&destlen_qual,&destlen_seq,(uint8_t*)data_ptr[i],len_arr[i],data_len[i],block_read_num[i]);
        else
            comtxt.COMTXT_decompress_unalign(pdest_qual[i],psrc_seq[i],&destlen_qual,&destlen_seq,(uint8_t*)data_ptr[i],len_arr[i],data_len[i],block_read_num[i],block_size[i]);
    }

    m_qualElPtr->m_rawBuf.updatesize(m_seqElPtr->m_rawBuf.size());
    m_seqElPtr->m_rawBuf.updatesize(m_seqElPtr->m_rawBuf.size());
    checkData();
    return 0;
}

