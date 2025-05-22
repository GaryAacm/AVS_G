#include "clr.h"
#include "simple_model.h"
#include "COMTXTCompress.h"
#include <cstring>
#include <map>
#include <algorithm>

#define MAX_SYMBOL 45
#define NUM_SYMBOL_QUAL 30
#define NUM_SYMBOL_SEQ 5
#define NUM_SYMBOL_COMBINE 45
#define XX_SYMBOL (MAX_SYMBOL-NUM_SYMBOL_QUAL)
using namespace std;


int COMTXTCompress::COMTXT_compress(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int cols = len_arr[0];
    int sourceLen = rows*cols;
    char *quality_block = new char[sourceLen];
    memcpy(quality_block,source,sizeof(unsigned char)*sourceLen);//将source存储在qul_block里面
    char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);//碱基存储在base里面


//============================建立编解码映射表=============================//

    map<pair<int,int>,int> mapp;
    map<pair<int,int>,int> mapp2;
    for(int i=0;i<=4;i++)
    {
        int k=0;
        for(int j=31;j<=38;j++)
        {
            mapp[make_pair(i,j)]=k+i*8;
            k++;
        }
    }
//----------------------------数据处理-------------------------------------//

    int model_num0 = NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*16;
    int model_num =5 * 5 * 5 * 5 + model_num0 * 4;//总的全部模型

    char table[256];
    memset(table, 0, 256);
    table['A'] = 0;
    table['T'] = 1;
    table['C'] = 2;
    table['G'] = 3;
    table['N'] = 4;

    char min_1 = '~';
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        if(quality_block[tmp]<min_1)
            min_1 = quality_block[tmp];
    }

    //-----------------------encode mean value of qual------------------------//数值编码

    int nums = 200;
    int count[nums + 1] = {0,}; //统计各个区间内的质量分数频率
    uint8_t* mean = (uint8_t*)calloc(rows, sizeof(uint8_t));//rows个序列，每个就是代表一行的序列

    //FILE *ss=fopen("compress_seq_use_to_compress.txt","w");


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
            gap[gap_index] = (double)i / 2;
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


//---------------------------------------small quality-not combine-------------------------------------------//
    int T = 6;
    SIMPLE_MODEL<NUM_SYMBOL_QUAL + 1> *model_qual;
    model_qual = new SIMPLE_MODEL<NUM_SYMBOL_QUAL + 1>[model_num + 1];
    RangeCoder rc3;
    rc3.output(dest);
    rc3.StartEncode();
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

            quality_block[offset] = quality_block[offset] - min_1;//量化

            if(quality_block[offset]>=31&&quality_block[offset]<=38) continue;


            char cur_base = base_block[offset];//对应的碱基
            int J0 = table[cur_base];//当前碱基对应的数字

            
            if (j == 0)//列，就是处理前两个
            {
                model_idx = 1;//给第一个进行分配
            }
            else if (j == 1)//第二个用前面一个来作为上下文编码
            {
                char pre_base_1 = base_block[offset - 1];//前面一个
                int Q1 = quality_block[offset - 1];//前一个质量分数

                model_idx = Q1 * 20 + 5;//分配模型，5是当前面对J0的处理，Q1已经归一化，保证是0到38不超出范围
            }
            else if (j == 2)
            {
                int Q1 = quality_block[offset - 1];
                int Q2 = quality_block[offset - 2];
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 +  5 + NUM_SYMBOL_QUAL * 20;//后面加上的是前面模型的偏置
            }
            else//用前三个
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基


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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 5 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }
            
            model_idx = model_idx + model_num0 * mean[i];
            if (quality_block[offset]>T)
            {
                model_qual[model_idx].encodeSymbol(&rc3, quality_block[offset]-T);
            }
            else
            {
                model_qual[model_idx].encodeSymbol(&rc3, 0);
                model_qual[model_num].encodeSymbol(&rc3, quality_block[offset]);
            }
        }
    }
    rc3.FinishEncode();
    int sz3=rc3.size_out();
    dest += sz3;
    *destLen += sz3;
    //===============================combine compress=========================//

    
    SIMPLE_MODEL<NUM_SYMBOL_COMBINE> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE>[model_num + 1];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

    //=====================================2 Process==============================//
    int i = 0, offset, model_idx_combine,model_idx_seq,model_idx_qual;
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
           
            int set_start = 5;
            //====================对seq=================//
            if(j==0)
            {
                model_idx_seq=0;
            }
            else if(j==1)
            {
                char pre_base_1 = base_block[offset - 1];
                int J1 = table[pre_base_1];
                model_idx_seq = J1 + 1;
            }
            else if(j==2)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基
                model_idx_seq = (J1 + 1) * (J2+1);

            }
            else if(j==3)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                char pre_base_3 = base_block[offset - 3];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                model_idx_seq = (J1+1)*(J2+1)*(J3+1);
            }
            else if(j==4)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                char pre_base_3 = base_block[offset - 3];
                char pre_base_4 = base_block[offset - 4];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
                model_idx_seq = (J1+1)*(J2+1)*(J3+1)*(J4+1);
            }
            else
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                char pre_base_3 = base_block[offset - 3];
                char pre_base_4 = base_block[offset - 4];
                char pre_base_5 = base_block[offset - 5];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
                int J5 = table[pre_base_5];
                model_idx_seq = (J1+1)*(J2+1)*(J3+1)*(J4+1)*(J5+1);
            }

            //=================对qual=====================//

            if(j==0)
            {
                
                model_idx_qual = 0;
            }
            else if(j==1)
            {
                int Q1 = max(quality_block[offset - 1]-30,0);
                model_idx = Q1 * 8 + 1;
            }
            else if(j==2)
            {

                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 8 +  1 + NUM_SYMBOL_QUAL * 8;
            }
            else
            {
                
                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                int Q3 = max(quality_block[offset - 3]-30,0);//前三个质量分数

                int A = Q1;
                int B = max(Q2,Q3);
                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = quality_block[offset - 4]-30;//前四个
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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 8 + B * 3  * 8 + C * 8  + 1 + NUM_SYMBOL_QUAL * 8 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 8;

            }
            

            //=========================== Encode Process =========================//
            model_idx_combine = (model_idx_qual) * model_idx_seq;
            if((quality_block[offset]>=31&&quality_block[offset]<=38))
            {
                model_combine[model_idx_combine].encodeSymbol(&rc,mapp[make_pair(int(table[base_block[offset]]),int(quality_block[offset]))]);
            }
            else
            {
                int seq_now = table[base_block[offset]];
                switch (seq_now)
                {
                case 0:
                    model_combine[model_idx_combine].encodeSymbol(&rc,40);
                    break;
                case 1:
                    model_combine[model_idx_combine].encodeSymbol(&rc,41); 
                    break;
                case 2:
                    model_combine[model_idx_combine].encodeSymbol(&rc,42);
                    break;
                case 3:
                    model_combine[model_idx_combine].encodeSymbol(&rc,43);
                    break;
                case 4:
                    model_combine[model_idx_combine].encodeSymbol(&rc,44);
                    break;
                default:
                    break;
                }
            }
        }
    }
    rc.FinishEncode();

	//fclose(ss);
    dest[rc.size_out()] = min_1;
    dest[rc.size_out()+1] = sz3;
    *destLen += rc.size_out() + 1;
    delete[] model_combine;
    delete[] model_qual;
    delete[] quality_block;
	delete[] base_block;
    return 0;
}

int COMTXTCompress::COMTXT_compress_unalign(char *dest,uint64_t *destLen,const uint8_t *source, const uint8_t *seq_source,const uint32_t *len_arr,uint64_t read_num)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int sourceLen = *destLen;
    char *quality_block = new char[sourceLen];
    memcpy(quality_block,source,sizeof(unsigned char)*sourceLen);
    //FILE *f =fopen("compress_seq_use_to_compress","w");
    char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);
    int seq_size = 1<<20;
    int COMTXT_MASK = seq_size - 1;

    int model_num =seq_size + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL;

    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
	// char *out0 = new char[100000000];
    char min_1 = '~';
//=================================组合编码的参照========================//

    map<pair<int,int>,int> mapp;
    for(int i=0;i<=4;i++)
    {
        int k=0;
        for(int j=31;j<=38;j++)
        {
            mapp[make_pair(i,j)]=k+i*8;
            k++;
        }
    }
//=============================================================================//

   //FILE *ss=fopen("compress_seq_use_to_compress.txt","w");

    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        if(quality_block[tmp]<min_1)
            min_1 = quality_block[tmp];
    }//一样的归一化

    int getclose[5]={0,6,12,18,24};
    //未对齐的质量分数
    // dest[0]=min_1;
    // memcpy(&dest[0],&min_1,sizeof(unsigned char));

    //===============================1.2 Encoder definition=========================//
	int T = XX_SYMBOL - 1;
	SIMPLE_MODEL<NUM_SYMBOL_COMBINE+5> *model_combine;
	model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE+5>[model_num + 4];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

    //=====================================2 Process==============================//
    int cur=0,offset=0, model_idx_seq,model_idx_qual,model_idx_combine;
    int last = 0x007616c7 & COMTXT_MASK; //设置先有的碱基上下文
    int last_start=last;
    //int *test_qual= new int[sourceLen],*test_seq=new int[sourceLen];
    //long long numm=0;

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

            if(j==0)
            {
                model_idx_seq = 0x007616c7 & COMTXT_MASK;
            }
            else
            {
                if(j<10)
                {
                    last = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=j-1;tmp>=0;tmp--)
                    {
                        int idx = offset-tmp-1;
                        char pre_base = base_block[idx];
                        uint8_t J1 = table[pre_base];
                        last=((last<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq = last;

                }
                else if(j==10)
                {
                    //last_start = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=9;tmp>=0;tmp--)
                    {
                        int idx=offset-tmp-1;
                        char pre_base = base_block[idx];
                        uint8_t J1 = table[pre_base];
                        last_start=((last_start<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq=last_start;
                }
                else
                {
                    char pre_base = base_block[offset-1];
                    uint8_t J1 = table[pre_base];
                    last_start=((last_start<<2)+J1) & COMTXT_MASK; //移动找到前10个
                    model_idx_seq=last_start;
                }
            }

            //===========================下面对qual=========================//

            if(j==0)
            {
                model_idx_qual = 0;
            }
            else if(j==1)
            {
                int Q1 = quality_block[offset - 1];//前一个质量分数

                model_idx_qual = Q1 * 20 +  1;
            }
            else if(j==2)
            {
                int Q1 = quality_block[offset - 1];
                int Q2 = quality_block[offset - 2];

                model_idx_qual = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 + 1 + NUM_SYMBOL_QUAL * 20;//后面加上的是前面模型的偏置

            }
            else
            {

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

                    //---------------------------------------------------------//
                    int A = Q1;//前一个
                    int B = max(Q2, Q3);//牵两个和三个的最大值
                    int C = 0;
                    //-----------------对质量分数-----------------------------//
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx_qual = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 1 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }

            model_idx_combine = model_idx_qual + model_idx_seq;

            //model_combine[model_idx_combine].encodeSymbol(&rc,table[base_block[offset]]);
            //uint16_t tt = table[base_block[offset]];
            //fprintf(ss,"%d\n",tt);

            //test_seq[numm] =model_idx_seq;
            //test_qual[numm]=model_idx_qual;
            //numm++;



            if(quality_block[offset]>=31&&quality_block[offset]<=38)
            {
                model_combine[model_idx_combine].encodeSymbol(&rc,mapp[make_pair(int(table[base_block[offset]]),int(quality_block[offset]))]);
            }

            else
            {
                if(quality_block[offset]<=5)
                {
                    model_combine[model_idx_combine].encodeSymbol(&rc, 40);
                    model_combine[model_num].encodeSymbol(&rc, table[base_block[offset]]);
                    model_combine[model_num].encodeSymbol(&rc, quality_block[offset]);
                }
                else if(quality_block[offset]<=11)
                {
                    model_combine[model_idx_combine].encodeSymbol(&rc, 41);
                    model_combine[model_num+1].encodeSymbol(&rc, table[base_block[offset]]+getclose[1]);
                    model_combine[model_num+1].encodeSymbol(&rc, quality_block[offset]);
                }
                else if(quality_block[offset]<=17)
                {
                    model_combine[model_idx_combine].encodeSymbol(&rc, 42);
                    model_combine[model_num+2].encodeSymbol(&rc, table[base_block[offset]]+getclose[2]);
                    model_combine[model_num+2].encodeSymbol(&rc, quality_block[offset]);
                }
                else if(quality_block[offset]<=23)
                {
                    model_combine[model_idx_combine].encodeSymbol(&rc, 43);
                    model_combine[model_num+3].encodeSymbol(&rc, table[base_block[offset]]+getclose[3]);
                    model_combine[model_num+3].encodeSymbol(&rc, quality_block[offset]);
                }
                else
                {
                    model_combine[model_idx_combine].encodeSymbol(&rc, 44);
                    model_combine[model_num+4].encodeSymbol(&rc, table[base_block[offset]]+getclose[4]);
                    model_combine[model_num+4].encodeSymbol(&rc, quality_block[offset]);
                }

            }

        }
    }
    rc.FinishEncode();
    dest[rc.size_out()] = min_1;
    *destLen = rc.size_out() + 1;

    delete[] model_combine;
    delete[] quality_block;
    delete[] base_block;
    return 0;
}



int COMTXTCompress::COMTXT_decompress(uint8_t *dest_qual,uint8_t *dest_seq,uint64_t *destLen_qual,uint64_t *destLen_seq,const uint8_t *source, const uint32_t *len_arr,int inlen,uint64_t read_num)
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
    char *in_buf_mean = new char[mean_sz];//存储质量分数均值序列
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

    char *decoded_matrix_qual = new char[sourceLen];//存储解码序列
    char *decoded_matrix_seq = new char[sourceLen];
    int seq_size=1<<20;
    int COMTXT_MASK = seq_size-1;
    int last = 0x007616c7 & COMTXT_MASK; //设置先有的碱基上下文
    int last_start=last;

    char *in_buf1 = new char[inlen-1];
    memcpy(in_buf1,source,sizeof(unsigned char)*(inlen-1));//把编码的数据存入，用于后面解码


    int model_num0 = NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL;
    int model_num = seq_size +  model_num0 * 4;//总的全部模型
    char table[256];
    memset(table, 0, 256);
	table['A'] = 0;
	table['T'] = 1;
	table['C'] = 2;
	table['G'] = 3;
	table['N'] = 4;
    char min_1 = source[inlen-1];//同编码

//=============================decode ref===============================//

    int seq[]={0,0,0,0,0,0,0,0,
               1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,
               4,4,4,4,4,4,4,4};


    int qual[]={31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38};

//=========================================================================//

    //===============================1.2 Encoder definition=========================//

    int T = XX_SYMBOL - 1; //6，低质量分数
    SIMPLE_MODEL<NUM_SYMBOL_COMBINE + 5> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE + 5>[model_num + 4];
    RangeCoder rc;
    rc.input(in_buf1); //初始化解码源数据
    rc.StartDecode();

    //=====================================2 Process==============================//
    int i = 0, offset, model_idx_qual,model_idx_seq,model_idx_combine;

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
            offset = i*cols + j;  //对应的行列的位置

//===================================start decode=======================================//


//===================================decode combine===================================//
 	    if(j==0)
            {
                model_idx_seq = last;
            }
            else
            {
                if(j<10)
                {
                    last = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=j-1;tmp>=0;tmp--)
                    {
                        int idx = offset-tmp-1;
                        char pre_base = decoded_matrix_seq[idx];
                        uint8_t J1 = table[pre_base];
                        last=((last<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq = last;

                }
                else if(j==10)
                {
                    //last_start = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=9;tmp>=0;tmp--)
                    {
                        int idx=offset-tmp-1;
                        char pre_base = decoded_matrix_seq[idx];
                        uint8_t J1 = table[pre_base];
                        last_start=((last_start<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq=last_start;
                }
                else
                {
                    char pre_base = decoded_matrix_seq[offset-1];
                    uint8_t J1 = table[pre_base];
                    last_start=((last_start<<2)+J1) & COMTXT_MASK; //移动找到前10个
                    model_idx_seq=last_start;
                }
            }

            //==================== qual decode ===========================//

            if(j==0)
            {
                model_idx_qual = 0;
            }
            else if(j==1)
            {
                int Q1 = decoded_matrix_qual[offset - 1];//前一个质量分数

                model_idx_qual = Q1 * 20 +  1;
            }
            else if(j==2)
            {
                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];

                model_idx_qual = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 + 1 + NUM_SYMBOL_QUAL * 20;//后面加上的是前面模型的偏置

            }
            else
            {

                    int Q1 = decoded_matrix_qual[offset - 1];
                    int Q2 = decoded_matrix_qual[offset - 2];
                    int Q3 = decoded_matrix_qual[offset - 3];//前三个质量分数

                    int Q4 = 0;
                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = decoded_matrix_qual[offset - 4];//前四个
                    }

                    //---------------------------------------------------------//
                    int A = Q1;//前一个
                    int B = max(Q2, Q3);//牵两个和三个的最大值
                    int C = 0;
                    //-----------------对质量分数-----------------------------//
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx_qual = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 1 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }

        model_idx_qual = model_idx_qual + model_num0 * mean[i];

        model_idx_combine = model_idx_qual+model_idx_seq;

	    int tmp=model_combine[model_idx_combine].decodeSymbol(&rc);

	if(tmp==40)
	{
	    decoded_matrix_seq[offset] = model_combine[model_num].decodeSymbol(&rc);
	    decoded_matrix_qual[offset] = model_combine[model_num].decodeSymbol(&rc);

	}
    else if(tmp==41)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+1].decodeSymbol(&rc)-getclose[1];
	    decoded_matrix_qual[offset] = model_combine[model_num+1].decodeSymbol(&rc);
    }
    else if(tmp==42)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+2].decodeSymbol(&rc)-getclose[2];
	    decoded_matrix_qual[offset] = model_combine[model_num+2].decodeSymbol(&rc);
    }
    else if(tmp==43)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+3].decodeSymbol(&rc)-getclose[3];
	    decoded_matrix_qual[offset] = model_combine[model_num+3].decodeSymbol(&rc);
    }
    else if(tmp==44)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+4].decodeSymbol(&rc)-getclose[4];
	    decoded_matrix_qual[offset] = model_combine[model_num+4].decodeSymbol(&rc);
    }
    else //combine compress
    {
        decoded_matrix_seq[offset] = seq[tmp];
        decoded_matrix_qual[offset] = qual[tmp];
    }

        if(int(decoded_matrix_seq[offset])==0) decoded_matrix_seq[offset]='A';
        else if(int(decoded_matrix_seq[offset])==1) decoded_matrix_seq[offset]='T';
        else if(int(decoded_matrix_seq[offset])==2) decoded_matrix_seq[offset]='C';
        else if(int(decoded_matrix_seq[offset])==3) decoded_matrix_seq[offset]='G';
        else if(int(decoded_matrix_seq[offset])==4) decoded_matrix_seq[offset]='N';

            /*if (int(decoded_matrix_qual[offset]) == 0)
            {
                decoded_matrix_qual[offset] = model_combine[model_num].decodeSymbol(&rc); //质量分数低的用最后一个数来存储
            }
            else
            {
                decoded_matrix_qual[offset] = decoded_matrix_qual[offset] + T;//高质量分数
            }*/
        }
    }
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix_qual[tmp] = decoded_matrix_qual[tmp] + min_1;
    }
    //归一化复原

    rc.FinishDecode();
    memcpy(dest_qual,decoded_matrix_qual,sizeof(unsigned char)*sourceLen);
    *destLen_qual = sourceLen;

    memcpy(dest_seq,decoded_matrix_seq,sizeof(unsigned char)*sourceLen);
    *destLen_seq = sourceLen;
    delete[] model_combine;
    delete[] in_buf1;
    delete[] decoded_matrix_qual;
     delete[] decoded_matrix_seq;

    delete[] in_buf_mean;
    free(mean);
    return 0;
}

int COMTXTCompress::COMTXT_decompress_unalign(uint8_t *dest_qual, uint8_t *dest_seq, uint64_t *destLen_qual, uint64_t *destLen_seq, const uint8_t *source, const uint32_t *len_arr, int inlen, uint64_t read_num, uint64_t size_file)
{
    //===============================1.1 Parameter definition & Initialization=======//
    int rows = read_num;
    int sourceLen = size_file;

    char *decoded_matrix_qual = new char[sourceLen];
    char *decoded_matrix_seq = new char[sourceLen];

    char *in_buf1 = new char[inlen-1];
    memcpy(in_buf1,source,sizeof(unsigned char)*(inlen-1));

    int seq_size=1<<20;
    int COMTXT_MASK = seq_size-1;

    int model_num = seq_size + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL;

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
    SIMPLE_MODEL<NUM_SYMBOL_COMBINE + 5> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE + 5>[model_num + 4];

    int getclose[5]={0,6,12,18,24};

    int seq[]={0,0,0,0,0,0,0,0,
               1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,
               4,4,4,4,4,4,4,4};


    int qual[]={31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38,
                31,32,33,34,35,36,37,38};


    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    //=====================================2 Process==============================//

    int cur = 0,offset = 0,model_idx_qual,model_idx_seq,model_idx_combine;
    int last = 0x007616c7 & COMTXT_MASK;
    int last_start=last;

    //FILE *f=fopen("decompress_seq_idx","w");
    //FILE *ff=fopen("decompress_idx_qual","w");
    //FILE *fff=fopen("decompress_use_to_test_seq_qual","w");

    for (int i = 0; i< rows; i++)
    {
        cur += len_arr[i];
        for (int j = 0; j < len_arr[i]; j++)  //列数
        {
            if(i==0)
                offset = j;
            else
                offset = cur-len_arr[i] + j;


           //===========================decode combine=====================================//

            if(j==0)
            {
                model_idx_seq = 0x007616c7 & COMTXT_MASK;
            }
            else
            {
                if(j<10)
                {
                    last = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=j-1;tmp>=0;tmp--)
                    {
                        int idx = offset-tmp-1;
                        char pre_base = decoded_matrix_seq[idx];
                        uint8_t J1 = table[pre_base];
                        last=((last<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq = last;

                }
                else if(j==10)
                {
                    //last_start = 0x007616c7 & COMTXT_MASK;
                    for(int tmp=9;tmp>=0;tmp--)
                    {
                        int idx=offset-tmp-1;
                        char pre_base = decoded_matrix_seq[idx];
                        uint8_t J1 = table[pre_base];
                        last_start=((last_start<<2) + J1) & COMTXT_MASK;
                    }
                    model_idx_seq=last_start;
                }
                else
                {
                    char pre_base = decoded_matrix_seq[offset-1];
                    uint8_t J1 = table[pre_base];
                    last_start=((last_start<<2)+J1) & COMTXT_MASK; //移动找到前10个
                    model_idx_seq=last_start;
                }
            }

            //===========================qual=========================//

            if(j==0)
            {
                model_idx_qual = 0;
            }
            else if(j==1)
            {
                int Q1 = decoded_matrix_qual[offset - 1];//前一个质量分数

                model_idx_qual = Q1 * 20 +  1;
            }
            else if(j==2)
            {
                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];

                model_idx_qual = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 + 1 + NUM_SYMBOL_QUAL * 20;//后面加上的是前面模型的偏置

            }
            else
            {

                    int Q1 = decoded_matrix_qual[offset - 1];
                    int Q2 = decoded_matrix_qual[offset - 2];
                    int Q3 = decoded_matrix_qual[offset - 3];//前三个质量分数

                    int Q4 = 0;
                    if (j == 3)
                    {
                        Q4 = Q3;
                    }
                    else
                    {
                        int Q4 = decoded_matrix_qual[offset - 4];//前四个
                    }

                    //---------------------------------------------------------//
                    int A = Q1;//前一个
                    int B = max(Q2, Q3);//牵两个和三个的最大值
                    int C = 0;
                    //-----------------对质量分数-----------------------------//
                    if (Q2 == Q3)
                    {
                        C = 1;
                    }
                    if (Q3 == Q4)
                    {
                        C = 2;
                    }
                    model_idx_qual = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 1 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }

            model_idx_combine = model_idx_qual + model_idx_seq;


	int tmp=int(model_combine[model_idx_combine].decodeSymbol(&rc));

	if(tmp==40)
	{
	    decoded_matrix_seq[offset] = model_combine[model_num].decodeSymbol(&rc);
	    decoded_matrix_qual[offset] = model_combine[model_num].decodeSymbol(&rc);

	}
    else if(tmp==41)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+1].decodeSymbol(&rc)-getclose[1];
	    decoded_matrix_qual[offset] = model_combine[model_num+1].decodeSymbol(&rc);
    }
    else if(tmp==42)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+2].decodeSymbol(&rc)-getclose[2];
	    decoded_matrix_qual[offset] = model_combine[model_num+2].decodeSymbol(&rc);
    }
    else if(tmp==43)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+3].decodeSymbol(&rc)-getclose[3];
	    decoded_matrix_qual[offset] = model_combine[model_num+3].decodeSymbol(&rc);
    }
    else if(tmp==44)
    {
        decoded_matrix_seq[offset] = model_combine[model_num+4].decodeSymbol(&rc)-getclose[4];
	    decoded_matrix_qual[offset] = model_combine[model_num+4].decodeSymbol(&rc);
    }
    else //combine compress
    {
        decoded_matrix_seq[offset] = seq[tmp];
        decoded_matrix_qual[offset] = qual[tmp];
    }



            if(int(decoded_matrix_seq[offset])==0) decoded_matrix_seq[offset]='A';
            else if(int(decoded_matrix_seq[offset])==1) decoded_matrix_seq[offset]='T';
            else if(int(decoded_matrix_seq[offset])==2) decoded_matrix_seq[offset]='C';
            else if(int(decoded_matrix_seq[offset])==3) decoded_matrix_seq[offset]='G';
            else if(int(decoded_matrix_seq[offset])==4) decoded_matrix_seq[offset]='N';

            /*if (int(decoded_matrix_qual[offset]) == 0)
            {
                decoded_matrix_qual[offset] = model_combine[model_num].decodeSymbol(&rc); //质量分数低的用最后一个数来存储
            }
            else
            {
                decoded_matrix_qual[offset] = decoded_matrix_qual[offset] + T;//高质量分数
            }*/

        }

    }

    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix_qual[tmp] = decoded_matrix_qual[tmp] + min_1;
    }

    rc.FinishDecode();


    //fclose(f);

    //fclose(ff);
    //fclose(fff);

    memcpy(dest_qual,decoded_matrix_qual,sizeof(unsigned char)*sourceLen);
    *destLen_qual = sourceLen;

    memcpy(dest_seq,decoded_matrix_seq,sizeof(unsigned char)*sourceLen);
    *destLen_seq = sourceLen;
    delete[] model_combine;
    delete[] in_buf1;
    delete[] decoded_matrix_qual;
     delete[] decoded_matrix_seq;

    return 0;
}
