#include "clr.h"
#include "simple_model.h"
#include "AcoCompress.h"
#include <cstring>
#include <algorithm>

#define MAX_SYMBOL 45
#define NUM_SYMBOL 38
#define XX_SYMBOL (MAX_SYMBOL-NUM_SYMBOL)


using namespace std;
// This function is used to compress the aligned quality score file 
int AcoCompress::aco_compress(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int cols = len_arr[0];
    int sourceLen = rows*cols;
    char *quality_block = new char[sourceLen];
    memcpy(quality_block,source,sizeof(unsigned char)*sourceLen);//将source存储在qul_block里面
    char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);//碱基存储在base里面
    int model_num0 = NUM_SYMBOL*NUM_SYMBOL*NUM_SYMBOL * 16;
    int model_num = model_num0 * 4;//总的全部模型
    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
	// char *out0 = new char[100000000];
    char min_1 = '~';
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        if(quality_block[tmp]<min_1)
            min_1 = quality_block[tmp];
    }//找到最小的质量分数

    //-----------------------encode mean value of qual------------------------//数值编码

    int nums = 200;
    int count[nums + 1] = {0,}; //统计各个区间内的质量分数频率
    uint8_t* mean = (uint8_t*)calloc(rows, sizeof(uint8_t));//rows个序列，每个就是代表一行的序列


    for (uint32_t i = 0; i < rows; ++i)
    {
        double sum = 0.0;
        for (uint32_t j = 0; j < cols; ++j)//质量分数
        {
            sum += source[i*cols + j];//按照列来读取来计算平均质量分数
        }
        sum /= cols;//平均值

        mean[i] = sum - min_1;//归一化，相对化

        int index = (sum - min_1) * 2;
        if (index < nums)
        {
            count[index]++;//计算此概率分布出现次数
        }
        else
        {
            count[nums]++;
        }
    }


    //用于对质量分数实行分区压缩
    double gap[3] = {0,};
    int gap_index = 0;
    int gap_sum = 0;
    for (uint32_t i = 0; i < nums + 1; i++)
    {
        gap_sum += count[i];
        if (gap_sum > read_num / 4)//超出序列长度四分之一
        {
            gap[gap_index] = (double)i / 2;//以此来分离质量分数高和低的区间，除以2对应前面的归一化*2
            gap_index++;
            gap_sum -= read_num / 4;
            if (gap_index == 3)
            {
                break;
            }
        }
    }


//按照列，将质量分数进行区间划分
    for (uint32_t i = 0; i < rows; ++i)
    {
        if (mean[i] < gap[0])
            mean[i] = 0;
        else if (mean[i] < gap[1])
            mean[i] = 1;
        else if (mean[i] < gap[2])
            mean[i] = 2;
        else
            mean[i] = 3;
    }

    //--------------coding-----------------//
    char *psize = dest;//输出码流，记录位置顺序
    dest += 4;//预留四个为了后面封装
    
	SIMPLE_MODEL<4> model_mean;
    RangeCoder rc2;
    rc2.output(dest);
    rc2.StartEncode();
    for (uint32_t i = 0; i < rows; ++i)
    {
        model_mean.encodeSymbol(&rc2, mean[i]);//对列质量分数的大小进行编码
    }
    rc2.FinishEncode();
    int sz = rc2.size_out();
    dest += sz;//压缩后的大小
    *destLen = 4 + sz;

    //------------Encap
    uint64_t sz_64 = (uint64_t)sz;
    sz_64 |= (uint64_t)1<<(4*7);
    for (int i = 0; i < 4; i++)
    {
        psize[i] = (sz_64>>(4-i-1)*8) & 0xff;
    }


//------------------------------------------------------------------------------------------
    //===============================1.2 Encoder definition=========================//
	int T = XX_SYMBOL - 1;
	SIMPLE_MODEL<NUM_SYMBOL+1> *model_qual;
	model_qual = new SIMPLE_MODEL<NUM_SYMBOL+1>[model_num + 1];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

    //=====================================2 Process==============================//
    int i = 0, offset, model_idx;
    for (int j = 0; j < cols; j++)
    {
        for (int k = 0; k < rows; k++)//行
        {
            if (j % 2 == 0)
            {
                i = k;//偶数位置
            }
            else
            {
                i = rows - k - 1;//之字形
            }

            offset = i*cols + j;//这一行计算当前元素在一维数组中的位置

            quality_block[offset] = quality_block[offset] - min_1;//计算归一化,为了保证从0开始到38

            char cur_base = base_block[offset];//对应的碱基
            int J0 = table[cur_base];//当前碱基对应的数字

            if (J0 == 4)//N碱基
            {
                model_idx = 0;
            }
            else
            {
                if (j == 0)//列，就是处理前两个
                {
                    model_idx = J0 + 1;//给第一个进行分配
                }
                else if (j == 1)//第二个用前面一个来作为上下文编码
                {
                    char pre_base_1 = base_block[offset - 1];//前面一个
                    int J1 = table[pre_base_1];

                    int G1 = J0 * 5 + J1;//用碱基计算，将当前碱基编码乘以5（因为碱基种类有4种，加上一个用于表示未知的N，共需要5个值来表示），再加上前一个碱基的编码，从而生成一个唯一的标识符
                    int Q1 = quality_block[offset - 1];//前一个质量分数

                    model_idx = Q1 * 20 + G1 + 5;//分配模型，5是当前面对J0的处理，Q1已经归一化，保证是0到38不超出范围
                }
                else if (j == 2)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = quality_block[offset - 1];
                    int Q2 = quality_block[offset - 2];
                    model_idx = (Q1*NUM_SYMBOL + Q2) * 20 + G1 + 5 + NUM_SYMBOL * 20;//后面加上的是前面模型的偏置
                }
                else//用前三个
                {
                    char pre_base_1 = base_block[offset - 1];
                    char pre_base_2 = base_block[offset - 2];
                    int J1 = table[pre_base_1];
                    int J2 = table[pre_base_2];//前两个碱基

                    int G1 = J0 * 5 + J1;//用前一个
                    int G2 = J0 * 4 * 4 + J1 * 4 + J2;//用前两个，不用考虑N，因为用了就不用考虑N

                    int Q1 = quality_block[offset - 1];
                    int Q2 = quality_block[offset - 2];
                    int Q3 = quality_block[offset - 3];//前三个质量分数

                    int Q4 = 0;
                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = quality_block[offset - 4];//前四个
                    }

                    int A = Q1;//前一个
                    int B = max(Q2, Q3);//牵两个和三个的最大值
                    int C = 0;
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    //这里把前四个的异同也作为考虑对象
                    model_idx = A * NUM_SYMBOL * 3  * 20 + B * 3  * 20 + C * 20  + G1 + 5 + NUM_SYMBOL * 20 + NUM_SYMBOL*NUM_SYMBOL * 20;
                }
            }
            model_idx = model_idx + model_num0 * mean[i];
            if (quality_block[offset]>T)
            {
                model_qual[model_idx].encodeSymbol(&rc, quality_block[offset]-T);
            }
            else
            {
                model_qual[model_idx].encodeSymbol(&rc, 0);
                model_qual[model_num].encodeSymbol(&rc, quality_block[offset]);
            }
        }
    }
    rc.FinishEncode();
	// int sz1;
    // sz1 = rc.size_out()+1;
    // memcpy(dest,out0,sizeof(unsigned char)*sz1);
    dest[rc.size_out()] = min_1;
    *destLen += rc.size_out() + 1;
    delete[] model_qual;
    delete[] quality_block;
	delete[] base_block;
    // delete[] out0;

    return 0;
}


//---------------------------------------------------------------------------------------------------//
// This function is used to decompress the aligned quality score file 
int AcoCompress::aco_decompress(uint8_t *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,int inlen,uint64_t read_num)
{
    //----------------------------decoding mean
    uint8_t* mean = (uint8_t*)calloc(read_num, sizeof(uint8_t));

    //-------------------decode size
    char *pbuf = (char*)source;
    uint32_t info_len = 0;
    uint64_t info_val;

    info_val = 0;
    uint8_t tmp = 0;
    for (int i = 7; i >= 0; i--)
    {
        tmp = 1<<i;
        if(pbuf[0] & tmp)
        {
            info_len = 8-i;
            tmp = pbuf[0] & (tmp - 1);
            info_val = tmp << (info_len-1)*8;
            break;
        }
    }

    int idx = 1;
    uint64_t tv;
    for (int i = info_len-2; i >= 0; i--,idx++)
    {
        tv = (uint8_t)pbuf[idx];
        info_val += tv<<(i*8);
    }
    int mean_sz = info_val;
    //-------------
    // int mean_sz = ACO_DECODE_INT((unsigned char *)source);
    source += 4;
    inlen -= 4;
    char *in_buf_mean = new char[mean_sz];//存储质量分数序列
    memcpy(in_buf_mean,source,sizeof(unsigned char)*(mean_sz));
    
	SIMPLE_MODEL<4> model_mean;//照应压缩，解码出四种情况对应不同质量分数的大小情况
    RangeCoder rc2;
    rc2.input(in_buf_mean);
    rc2.StartDecode();
    for (int i = 0; i < read_num; i++)
    {
        mean[i] = model_mean.decodeSymbol(&rc2);
    }
    rc2.FinishDecode();
    source += mean_sz;
    inlen -= mean_sz;

    //===============================1.1 Parameter definition & Initialization=======//
    int cols = len_arr[0];
    int rows = read_num;
    int sourceLen = rows*cols;

    char *decoded_matrix = new char[sourceLen];//存储解码序列
    char *in_buf1 = new char[inlen-1];
    memcpy(in_buf1,source,sizeof(unsigned char)*(inlen-1));

	char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);//存储碱基

    int model_num0 = NUM_SYMBOL*NUM_SYMBOL*NUM_SYMBOL * 16;
    int model_num = model_num0 * 4;//*4是为了区别质量分数不同的四类碱基
    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
    char min_1 = source[inlen-1];//同编码



    //===============================1.2 Encoder definition=========================//

    int T = XX_SYMBOL - 1; //6，低质量分数
    SIMPLE_MODEL<NUM_SYMBOL + 1> *model_qual;
	model_qual = new SIMPLE_MODEL<NUM_SYMBOL + 1>[model_num + 1];
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    //=====================================2 Process==============================//
    int i = 0, offset, model_idx;
    for (int j = 0; j < cols; j++)//操作同编码，找到对应模型
    {
        for (int k = 0; k < rows; k++)
        {
            if (j % 2 == 0)
            {
                i = k;
            }
            else
            {
                i = rows - k - 1;
            }
            offset = i*cols + j;
            char cur_base = base_block[offset];
            int J0 = table[cur_base];
            if (J0 == 4)
            {
                model_idx = 0;
            }
            else
            {
                if (j == 0)
                {
                    model_idx = J0 + 1;
                }
                else if (j == 1)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = decoded_matrix[offset - 1];
                    model_idx = Q1 * 20 + G1 + 5;
                }
                else if (j == 2)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = decoded_matrix[offset - 1];
                    int Q2 = decoded_matrix[offset - 2];
                    model_idx = (Q1*NUM_SYMBOL + Q2) * 20 + G1 + 5 + NUM_SYMBOL * 20;
                }
                else
                {
                    char pre_base_1 = base_block[offset - 1];
                    char pre_base_2 = base_block[offset - 2];
                    int J1 = table[pre_base_1];
                    int J2 = table[pre_base_2];
                    int G1 = J0 * 5 + J1;
                    int G2 = J0 * 4 * 4 + J1 * 4 + J2;
                    int Q1 = decoded_matrix[offset - 1];
                    int Q2 = decoded_matrix[offset - 2];
                    int Q3 = decoded_matrix[offset - 3];
                    int Q4 = 0;
                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = decoded_matrix[offset - 4];
                    }
                    int A = Q1;
                    int B = max(Q2, Q3);
                    int C = 0;
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx = A * NUM_SYMBOL * 3  * 20 + B * 3  * 20 + C * 20  + G1 + 5 + NUM_SYMBOL * 20 + NUM_SYMBOL*NUM_SYMBOL * 20;
                }
            }
            model_idx = model_idx + model_num0 * mean[i];

            //用同样的模型编解码
            decoded_matrix[offset] = model_qual[model_idx].decodeSymbol(&rc);
            if (int(decoded_matrix[offset]) == 0)
            {
                decoded_matrix[offset] = model_qual[model_num].decodeSymbol(&rc);
            }
            else
            {
                decoded_matrix[offset] = decoded_matrix[offset] + T;//高质量分数
            }
        }
    }
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix[tmp] = decoded_matrix[tmp] + min_1;
    }
    
    rc.FinishDecode();
    memcpy(dest,decoded_matrix,sizeof(unsigned char)*sourceLen);
    *destLen = sourceLen;//跟新数据长度
    delete[] model_qual;
    delete[] in_buf1;
    delete[] base_block;
    delete[] decoded_matrix;
    
    delete[] in_buf_mean;
    free(mean);
    return 0;
}  


// This function is used to compress the unaligned quality score file 
int AcoCompress::aco_compress_unalign(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int sourceLen = *destLen;
    char *quality_block = new char[sourceLen];
    memcpy(quality_block,source,sizeof(unsigned char)*sourceLen);
    char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);
    int model_num = NUM_SYMBOL*NUM_SYMBOL*NUM_SYMBOL * 16;
    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
	// char *out0 = new char[100000000];
    char min_1 = '~';
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        if(quality_block[tmp]<min_1)
            min_1 = quality_block[tmp];
    }//一样的归一化


    //未对齐的质量分数
    // dest[0]=min_1;
    // memcpy(&dest[0],&min_1,sizeof(unsigned char));

    //===============================1.2 Encoder definition=========================//
	int T = XX_SYMBOL - 1;
	SIMPLE_MODEL<NUM_SYMBOL+1> *model_qual;
	model_qual = new SIMPLE_MODEL<NUM_SYMBOL+1>[model_num + 1];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

    //=====================================2 Process==============================//
    int cur=0,offset=0, model_idx;
    for (int i = 0; i< rows; i++)
    {
        cur += len_arr[i];//碱基长度
        for (int j = 0; j < len_arr[i]; j++)//直接用每行的碱基长度去对应
        {
            if(i==0)
                offset = j;
            else
                offset = cur-len_arr[i] + j;//类似的计算偏移量

            quality_block[offset] = quality_block[offset] - min_1;
            char cur_base = base_block[offset];
            int J0 = table[cur_base];
            if (J0 == 4)//类似的N碱基
            {
                model_idx = 0;
            }
            else
            {
                if (j == 0)
                {
                    model_idx = J0 + 1;
                }
                else if (j == 1)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = quality_block[offset - 1];
                    model_idx = Q1 * 20 + G1 + 5;
                }
                else if (j == 2)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = quality_block[offset - 1];
                    int Q2 = quality_block[offset - 2];
                    model_idx = (Q1*NUM_SYMBOL + Q2) * 20 + G1 + 5 + NUM_SYMBOL * 20;
                }
                else
                {
                    char pre_base_1 = base_block[offset - 1];
                    char pre_base_2 = base_block[offset - 2];
                    int J1 = table[pre_base_1];
                    int J2 = table[pre_base_2];
                    int G1 = J0 * 5 + J1;
                    int G2 = J0 * 4 * 4 + J1 * 4 + J2;
                    int Q1 = quality_block[offset - 1];
                    int Q2 = quality_block[offset - 2];
                    int Q3 = quality_block[offset - 3];
                    int Q4 = 0;

                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = quality_block[offset - 4];
                    }
                    int A = Q1;
                    int B = max(Q2, Q3);
                    int C = 0;
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx = A * NUM_SYMBOL * 3  * 20 + B * 3  * 20 + C * 20 + G1 + 5 + NUM_SYMBOL * 20 + NUM_SYMBOL*NUM_SYMBOL * 20;
                }
            }
            if (quality_block[offset]>T)
            {
                model_qual[model_idx].encodeSymbol(&rc, quality_block[offset]-T);
            }
            else
            {
                model_qual[model_idx].encodeSymbol(&rc, 0);
                model_qual[model_num].encodeSymbol(&rc, quality_block[offset]);
            }
        }
    }
    rc.FinishEncode();
	// int sz1;
    // sz1 = rc.size_out()+1;
    // memcpy(dest,out0,sizeof(unsigned char)*sz1);
    dest[rc.size_out()] = min_1;
    *destLen = rc.size_out() + 1;
    delete[] model_qual;
    delete[] quality_block;
	delete[] base_block;
    // delete[] out0; 
    return 0;
}
//与对齐不同，用的每行碱基的序列去对应


// This function is used to decompress the unaligned quality score file 
int AcoCompress::aco_decompress_unalign(uint8_t *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,int inlen,uint64_t read_num,uint64_t size_file)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int sourceLen = size_file;
    char *decoded_matrix = new char[sourceLen];
    char *in_buf1 = new char[inlen-1];
    memcpy(in_buf1,source,sizeof(unsigned char)*(inlen-1));
	char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);
    int model_num = NUM_SYMBOL*NUM_SYMBOL*NUM_SYMBOL * 16;
    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
    char min_1 = source[inlen-1];
    //===============================1.2 Encoder definition=========================//
    int T = XX_SYMBOL - 1;
    SIMPLE_MODEL<NUM_SYMBOL + 1> *model_qual;
	model_qual = new SIMPLE_MODEL<NUM_SYMBOL + 1>[model_num + 1];
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    //=====================================2 Process==============================//
    int cur=0,offset=0, model_idx;
    for (int i = 0; i< rows; i++)
    {
        cur += len_arr[i];
        for (int j = 0; j < len_arr[i]; j++)//解码部分也只是序列长度处理的不同
        {
            if(i==0)
                offset = j;
            else
                offset = cur-len_arr[i] + j;

            char cur_base = base_block[offset];
            int J0 = table[cur_base];
            if (J0 == 4)
            {
                model_idx = 0;
            }
            else
            {
                if (j == 0)
                {
                    model_idx = J0 + 1;
                }
                else if (j == 1)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = decoded_matrix[offset - 1];
                    model_idx = Q1 * 20 + G1 + 5;
                }
                else if (j == 2)
                {
                    char pre_base_1 = base_block[offset - 1];
                    int J1 = table[pre_base_1];
                    int G1 = J0 * 5 + J1;
                    int Q1 = decoded_matrix[offset - 1];
                    int Q2 = decoded_matrix[offset - 2];
                    model_idx = (Q1*NUM_SYMBOL + Q2) * 20 + G1 + 5 +  * 20;
                }
                else
                {
                    char pre_base_1 = base_block[offset - 1];
                    char pre_base_2 = base_block[offset - 2];
                    int J1 = table[pre_base_1];
                    int J2 = table[pre_base_2];
                    int G1 = J0 * 5 + J1;
                    int G2 = J0 * 4 * 4 + J1 * 4 + J2;
                    int Q1 = decoded_matrix[offset - 1];
                    int Q2 = decoded_matrix[offset - 2];
                    int Q3 = decoded_matrix[offset - 3];
                    int Q4 = 0;
                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = decoded_matrix[offset - 4];
                    }
                    int A = Q1;
                    int B = max(Q2, Q3);
                    int C = 0;
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx = A * NUM_SYMBOL * 3  * 20 + B * 3  * 20 + C * 20 + G1 + 5 + NUM_SYMBOL * 20 + NUM_SYMBOL*NUM_SYMBOL * 20;
                }
            }
            decoded_matrix[offset] = model_qual[model_idx].decodeSymbol(&rc);
            if (int(decoded_matrix[offset]) == 0)
            {
                decoded_matrix[offset] = model_qual[model_num].decodeSymbol(&rc);
            }
            else
            {
                decoded_matrix[offset] = decoded_matrix[offset] + T;
            }
            // decoded_matrix[offset] = decoded_matrix[offset] + '!';
        }
    }
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix[tmp] = decoded_matrix[tmp] + min_1;
    }
    rc.FinishDecode();
    memcpy(dest,decoded_matrix,sizeof(unsigned char)*sourceLen);
    *destLen = sourceLen;
    delete[] model_qual;
    delete[] in_buf1;
    delete[] base_block;
    delete[] decoded_matrix;
    return 0;
}
