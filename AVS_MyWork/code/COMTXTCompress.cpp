#include "clr.h"
#include "simple_model.h"
#include "COMTXTCompress.h"
#include <cstring>
#include <map>
#include <algorithm>

#define MAX_SYMBOL 45
#define NUM_SYMBOL_QUAL 38
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
    int model_num =5 * 5 * 5 * 5 *5 * 8 * 8 * 8 *2 + model_num0 * 4;//总的全部模型
    int pianyi = 5 * 5 * 5 * 5 *5 * 8 * 8 * 8 *2 ;
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
    uint8_t* mean = (uint8_t*)calloc(rows, sizeof(uint8_t));

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
    char *psize = dest;
    dest += 4;

	SIMPLE_MODEL<4> model_mean;
    RangeCoder rc2;
    rc2.output(dest);
    rc2.StartEncode();
    for (uint32_t i = 0; i < rows; ++i)
    {
        model_mean.encodeSymbol(&rc2, mean[i]);
    }
    rc2.FinishEncode();
    int sz = rc2.size_out();
    dest += sz;
    *destLen = 4 + sz;

    //------------Encap
    uint64_t sz_64 = (uint64_t)sz;
    sz_64 |= (uint64_t)1<<(4*7);
    for (int i = 0; i < 4; i++)
    {
        psize[i] = (sz_64>>(4-i-1)*8) & 0xff;
    }


//---------------------------------------small quality-not combine-------------------------------------------//
    //===============================combine compress=========================//

    SIMPLE_MODEL<NUM_SYMBOL_COMBINE> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE>[model_num + 1];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

   //=====================================2 Process==============================//
    int model_idx_combine,model_idx_seq,model_idx_qual,model_idx,i=0,offset;
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
           
            int set_start = 5;
            //====================对seq=================//
            if(j==0)
            {
                model_idx_seq=5;
            }
            else if(j==1)
            {
                char pre_base_1 = base_block[offset - 1];
                int J1 = table[pre_base_1];
                model_idx_seq = J1+1+5;
            }
            else if(j==2)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基                
                model_idx_seq = (J1+1)*5+J2 + 1 + 5;

            }
            else if(j==3)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                char pre_base_3 = base_block[offset - 3];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
             
                model_idx_seq = (J1+1)*5*5+(J2+1)*5+J3+1+5;
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
                //model_idx_seq = (J1+1)*(J2+1)*(J3+1)*(J4+1)+5;
                model_idx_seq=(J1+1)*5*5*5+(J2+1)*5*5+(J3+1)*5+J4+1+5;
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
               
                model_idx_seq=(J1+1)*5*5*5*5+(J2+1)*5*5*5+(J3+1)*5*5+J4*5+J5+1+5;
            }

            //=================对qual=====================//

            if(j==0)
            {
                
                model_idx_qual = 1;
            }
            else if(j==1)
            {
                int Q1 = max(quality_block[offset - 1]-30,0);
                //model_idx_qual = Q1 + 1;
                model_idx_qual = (Q1+1) * 8 + 1;
            }
            else if(j==2)
            {

                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                //model_idx_qual = (Q1+1)*(Q2+1);
                model_idx_qual = ((Q1+1)*8 + Q2+1) * 8 +  1 + 8 * 8;
            }
            else
            {
                
                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                int Q3 = max(quality_block[offset - 3]-30,0);//前三个质量分数

                int A = Q1+1;
                int B = max(Q2,Q3)+1;
                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = max(quality_block[offset - 4]-30,0);//前四个
                }

               int C=0;
                if (Q2 == Q3)
                {
                    C = 1;
                }
                if (Q3 == Q4)
                {
                    C = 2;
                }
                
                model_idx_qual = A * 8 * 8 + B * 8 + C * 8  + 1 + 8 * 8 + 8 *8 * 8;

            }
            
	    //================================small quality  combine==================================//
	    if (j == 0)
            {
                model_idx = 1;
            }
            else if (j == 1)
            {
                char pre_base_1 = base_block[offset - 1];
                int Q1 = quality_block[offset - 1];
                model_idx = Q1 * 20 + 5;
            }
            else if (j == 2)
            {
                int Q1 = quality_block[offset - 1];
                int Q2 = quality_block[offset - 2];
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 +  5 + NUM_SYMBOL_QUAL * 20;
            }
            else
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];


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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 5 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }
            
            //=========================== Encode Process =========================//
            model_idx_combine = model_idx_qual * model_idx_seq;
            model_idx = model_idx + model_num0 * mean[i];
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
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 1:
                    model_combine[model_idx_combine].encodeSymbol(&rc,41);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 2:
                    model_combine[model_idx_combine].encodeSymbol(&rc,42);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 3:
                    model_combine[model_idx_combine].encodeSymbol(&rc,43);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 4:
                    model_combine[model_idx_combine].encodeSymbol(&rc,44);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                default:
                    break;
                }
            }
        }
    }
    rc.FinishEncode();

    dest[rc.size_out()] = min_1;
    *destLen += rc.size_out() + 1;
    delete[] model_combine;
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
    char *base_block = new char[sourceLen];
    memcpy(base_block,seq_source,sizeof(unsigned char)*sourceLen);

    int model_num = NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*16+5*5*5*5*5*8*8*8*2;
    int pianyi = 5*5*5*5*5*8*8*8*2;

    char table[256];
    memset(table, 0, 256);
table['A'] = 0;
table['T'] = 1;
table['C'] = 2;
table['G'] = 3;
table['N'] = 4;
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

    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        if(quality_block[tmp]<min_1)
            min_1 = quality_block[tmp];
    }

 
    //===============================1.2 Encoder definition=========================//
	int T = XX_SYMBOL - 1;
	SIMPLE_MODEL<NUM_SYMBOL_COMBINE> *model_combine;
	model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE>[model_num + 1];
    RangeCoder rc;
    rc.output(dest);
    rc.StartEncode();

    //=====================================2 Process==============================//
    int model_idx_combine,model_idx_seq,model_idx_qual,model_idx,i=0,offset,cur=0;
    for (int i = 0; i < rows; i++)
    {
        cur += len_arr[i];
        for (int j = 0; j < len_arr[i]; j++)
        {
            if(i==0)
                offset = j;
            else
                offset = cur-len_arr[i] + j;
            quality_block[offset] = quality_block[offset] - min_1;//量化
           
            int set_start = 5;
            //====================对seq=================//
            if(j==0)
            {
                model_idx_seq=5;
            }
            else if(j==1)
            {
                char pre_base_1 = base_block[offset - 1];
                int J1 = table[pre_base_1];
                model_idx_seq = J1+1+5;
            }
            else if(j==2)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基                
                model_idx_seq = (J1+1)*5+J2 + 1 + 5;

            }
            else if(j==3)
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                char pre_base_3 = base_block[offset - 3];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
             
                model_idx_seq = (J1+1)*5*5+(J2+1)*5+J3+1+5;
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
       
                model_idx_seq=(J1+1)*5*5*5+(J2+1)*5*5+(J3+1)*5+J4+1+5;
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
               
                model_idx_seq=(J1+1)*5*5*5*5+(J2+1)*5*5*5+(J3+1)*5*5+J4*5+J5+1+5;
            }

            //=================对qual=====================//

            if(j==0)
            {
                
                model_idx_qual = 1;
            }
            else if(j==1)
            {
                int Q1 = max(quality_block[offset - 1]-30,0);
                //model_idx_qual = Q1 + 1;
                model_idx_qual = (Q1+1) * 8 + 1;
            }
            else if(j==2)
            {

                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                //model_idx_qual = (Q1+1)*(Q2+1);
                model_idx_qual = ((Q1+1)*8 + Q2+1) * 8 +  1 + 8 * 8;
            }
            else
            {
                
                int Q1 = max(quality_block[offset - 1]-30,0);
                int Q2 = max(quality_block[offset - 2]-30,0);
                int Q3 = max(quality_block[offset - 3]-30,0);//前三个质量分数

                int A = Q1+1;
                int B = max(Q2,Q3)+1;
                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = max(quality_block[offset - 4]-30,0);//前四个
                }

               int C=0;
                if (Q2 == Q3)
                {
                    C = 1;
                }
                if (Q3 == Q4)
                {
                    C = 2;
                }
                
                model_idx_qual = A * 8 * 8 + B * 8 + C * 8  + 1 + 8 * 8 + 8 *8 * 8;

            }
            
	    //================================small combine==================================//
	    if (j == 0)
            {
                model_idx = 1;
            }
            else if (j == 1)
            {
                char pre_base_1 = base_block[offset - 1];
                int Q1 = quality_block[offset - 1];
                model_idx = Q1 * 20 + 5;
            }
            else if (j == 2)
            {
                int Q1 = quality_block[offset - 1];
                int Q2 = quality_block[offset - 2];
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 +  5 + NUM_SYMBOL_QUAL * 20;
            }
            else
            {
                char pre_base_1 = base_block[offset - 1];
                char pre_base_2 = base_block[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];


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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 5 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }
            
            //=========================== Encode Process =========================//
            model_idx_combine = model_idx_qual * model_idx_seq;
            model_idx = model_idx;
            
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
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 1:
                    model_combine[model_idx_combine].encodeSymbol(&rc,41);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 2:
                    model_combine[model_idx_combine].encodeSymbol(&rc,42);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 3:
                    model_combine[model_idx_combine].encodeSymbol(&rc,43);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                case 4:
                    model_combine[model_idx_combine].encodeSymbol(&rc,44);
                    model_combine[model_idx+pianyi].encodeSymbol(&rc,quality_block[offset]);
                    break;
                default:
                    break;
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
    source += 4;
    inlen -= 4;
    char *in_buf_mean = new char[mean_sz];//存储质量分数均值序列
    memcpy(in_buf_mean,source,sizeof(unsigned char)*(mean_sz));

    SIMPLE_MODEL<4> model_mean;
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

    char *in_buf1 = new char[inlen-1];
    memcpy(in_buf1,source,sizeof(unsigned char)*(inlen-1));//把编码的数据存入，用于后面解码


    int model_num0 = NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*16;
    int model_num = 5*5*5*5*5*8*8*8*2 +  model_num0 * 4;//总的全部模型
    int pianyi = 5*5*5*5*5*8*8*8*2 ;
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

    int T = XX_SYMBOL - 1; 
    SIMPLE_MODEL<NUM_SYMBOL_COMBINE> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE>[model_num + 1];
    RangeCoder rc;
    rc.input(in_buf1);
    rc.StartDecode();

    //=====================================2 Process==============================//
    int i = 0, offset, model_idx_qual,model_idx_seq,model_idx_combine,model_idx;

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

//================================== decode combine=======================================//
	   if(j==0)
            {
                model_idx_seq=5;
            }
            else if(j==1)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                int J1 = table[pre_base_1];
                model_idx_seq = J1+1+5;
            }
            else if(j==2)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基                
                model_idx_seq = (J1+1)*5+J2 + 1 + 5;

            }
            else if(j==3)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
             
                model_idx_seq = (J1+1)*5*5+(J2+1)*5+J3+1+5;
            }
            else if(j==4)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                char pre_base_4 = decoded_matrix_seq[offset - 4];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
           
                model_idx_seq=(J1+1)*5*5*5+(J2+1)*5*5+(J3+1)*5+J4+1+5;
            }
            else
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                char pre_base_4 = decoded_matrix_seq[offset - 4];
                char pre_base_5 = decoded_matrix_seq[offset - 5];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
                int J5 = table[pre_base_5];
               
                model_idx_seq=(J1+1)*5*5*5*5+(J2+1)*5*5*5+(J3+1)*5*5+J4*5+J5+1+5;
            }

            //=================对qual=====================//

            if(j==0)
            {
                
                model_idx_qual = 1;
            }
            else if(j==1)
            {
                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                model_idx_qual = (Q1+1) * 8 + 1;
            }
            else if(j==2)
            {

                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                int Q2 = max(decoded_matrix_qual[offset - 2]-30,0);
                //model_idx_qual = (Q1+1)*(Q2+1);
                model_idx_qual = ((Q1+1)*8 + Q2+1) * 8 +  1 + 8 * 8;
            }
            else
            {
                
                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                int Q2 = max(decoded_matrix_qual[offset - 2]-30,0);
                int Q3 = max(decoded_matrix_qual[offset - 3]-30,0);//前三个质量分数

                int A = Q1+1;
                int B = max(Q2,Q3)+1;
                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = max(decoded_matrix_qual[offset - 4]-30,0);//前四个
                }

               int C=0;
                if (Q2 == Q3)
                {
                    C = 1;
                }
                if (Q3 == Q4)
                {
                    C = 2;
                }
                
                model_idx_qual = A * 8 * 8 + B * 8 + C * 8  + 1 + 8 * 8 + 8 *8 * 8;

            }
            
            //================================small quality  combine==================================//
	    if (j == 0)
            {
                model_idx = 1;
            }
            else if (j == 1)
            {
                char pre_base_1 = decoded_matrix_qual[offset - 1];
                int Q1 = decoded_matrix_qual[offset - 1];
                model_idx = Q1 * 20 + 5;
            }
            else if (j == 2)
            {
                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 +  5 + NUM_SYMBOL_QUAL * 20;
            }
            else
            {
                char pre_base_1 = decoded_matrix_qual[offset - 1];
                char pre_base_2 = decoded_matrix_qual[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];


                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];
                int Q3 = decoded_matrix_qual[offset - 3];

                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = decoded_matrix_qual[offset - 4];
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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 5 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }
 	     model_idx_combine = model_idx_qual * model_idx_seq;
            model_idx = model_idx + model_num0 * mean[i];
            int tmp=model_combine[model_idx_combine].decodeSymbol(&rc);
            if(tmp==40)
	{
	    decoded_matrix_seq[offset] = 'A';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);

	}
    else if(tmp==41)
    {
        decoded_matrix_seq[offset] = 'T';
	decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==42)
    {
        decoded_matrix_seq[offset] = 'C';
	decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==43)
    {
        decoded_matrix_seq[offset] = 'G';
	decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==44)
    {
        decoded_matrix_seq[offset] = 'N';
	decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else //combine compress
    {
        decoded_matrix_seq[offset] = seq[tmp];
        decoded_matrix_qual[offset] = qual[tmp];
        if(int(decoded_matrix_seq[offset])==0) decoded_matrix_seq[offset]='A';
        else if(int(decoded_matrix_seq[offset])==1) decoded_matrix_seq[offset]='T';
        else if(int(decoded_matrix_seq[offset])==2) decoded_matrix_seq[offset]='C';
        else if(int(decoded_matrix_seq[offset])==3) decoded_matrix_seq[offset]='G';
        else if(int(decoded_matrix_seq[offset])==4) decoded_matrix_seq[offset]='N';
    }
    
   }
  }
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix_qual[tmp] = decoded_matrix_qual[tmp] + min_1;
    }
   

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


    int model_num = 5*5*5*5*5*8*8*8*2 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL*16;
    int pianyi=5*5*5*5*5*8*8*8*2;

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
    SIMPLE_MODEL<NUM_SYMBOL_COMBINE> *model_combine;
    model_combine = new SIMPLE_MODEL<NUM_SYMBOL_COMBINE>[model_num + 1];


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

    int cur = 0,offset = 0,model_idx_qual,model_idx_seq,model_idx_combine,model_idx;
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
                model_idx_seq=5;
            }
            else if(j==1)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                int J1 = table[pre_base_1];
                model_idx_seq = J1+1+5;
            }
            else if(j==2)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];//前两个碱基                
                model_idx_seq = (J1+1)*5+J2 + 1 + 5;

            }
            else if(j==3)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
             
                model_idx_seq = (J1+1)*5*5+(J2+1)*5+J3+1+5;
            }
            else if(j==4)
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                char pre_base_4 = decoded_matrix_seq[offset - 4];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
           
                model_idx_seq=(J1+1)*5*5*5+(J2+1)*5*5+(J3+1)*5+J4+1+5;
            }
            else
            {
                char pre_base_1 = decoded_matrix_seq[offset - 1];
                char pre_base_2 = decoded_matrix_seq[offset - 2];
                char pre_base_3 = decoded_matrix_seq[offset - 3];
                char pre_base_4 = decoded_matrix_seq[offset - 4];
                char pre_base_5 = decoded_matrix_seq[offset - 5];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];
                int J3 = table[pre_base_3];
                int J4 = table[pre_base_4];
                int J5 = table[pre_base_5];
               
                model_idx_seq=(J1+1)*5*5*5*5+(J2+1)*5*5*5+(J3+1)*5*5+J4*5+J5+1+5;
            }

            //=================对qual=====================//

            if(j==0)
            {
                
                model_idx_qual = 1;
            }
            else if(j==1)
            {
                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                model_idx_qual = (Q1+1) * 8 + 1;
            }
            else if(j==2)
            {

                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                int Q2 = max(decoded_matrix_qual[offset - 2]-30,0);
                //model_idx_qual = (Q1+1)*(Q2+1);
                model_idx_qual = ((Q1+1)*8 + Q2+1) * 8 +  1 + 8 * 8;
            }
            else
            {
                
                int Q1 = max(decoded_matrix_qual[offset - 1]-30,0);
                int Q2 = max(decoded_matrix_qual[offset - 2]-30,0);
                int Q3 = max(decoded_matrix_qual[offset - 3]-30,0);//前三个质量分数

                int A = Q1+1;
                int B = max(Q2,Q3)+1;
                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = max(decoded_matrix_qual[offset - 4]-30,0);//前四个
                }

               int C=0;
                if (Q2 == Q3)
                {
                    C = 1;
                }
                if (Q3 == Q4)
                {
                    C = 2;
                }
                
                model_idx_qual = A * 8 * 8 + B * 8 + C * 8  + 1 + 8 * 8 + 8 *8 * 8;

            }
            
            //================================small quality  combine==================================//
	    if (j == 0)
            {
                model_idx = 1;
            }
            else if (j == 1)
            {
                char pre_base_1 = decoded_matrix_qual[offset - 1];
                int Q1 = decoded_matrix_qual[offset - 1];
                model_idx = Q1 * 20 + 5;
            }
            else if (j == 2)
            {
                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];
                model_idx = (Q1*NUM_SYMBOL_QUAL + Q2) * 20 +  5 + NUM_SYMBOL_QUAL * 20;
            }
            else
            {
                char pre_base_1 = decoded_matrix_qual[offset - 1];
                char pre_base_2 = decoded_matrix_qual[offset - 2];
                int J1 = table[pre_base_1];
                int J2 = table[pre_base_2];


                int Q1 = decoded_matrix_qual[offset - 1];
                int Q2 = decoded_matrix_qual[offset - 2];
                int Q3 = decoded_matrix_qual[offset - 3];

                int Q4 = 0;
                if (j == 3)
                {
                    Q4 = Q3;
                }
                else
                {
                    int Q4 = decoded_matrix_qual[offset - 4];
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
                model_idx = A * NUM_SYMBOL_QUAL * 3  * 20 + B * 3  * 20 + C * 20  + 5 + NUM_SYMBOL_QUAL * 20 + NUM_SYMBOL_QUAL*NUM_SYMBOL_QUAL * 20;
            }
 	     model_idx_combine = model_idx_qual * model_idx_seq;
            
            int tmp=model_combine[model_idx_combine].decodeSymbol(&rc);
            if(tmp==40)
	{
	    decoded_matrix_seq[offset] = 'A';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);

	}
    else if(tmp==41)
    {
        decoded_matrix_seq[offset] = 'T';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==42)
    {
        decoded_matrix_seq[offset] = 'C';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==43)
    {
        decoded_matrix_seq[offset] = 'G';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else if(tmp==44)
    {
        decoded_matrix_seq[offset] = 'N';
	    decoded_matrix_qual[offset] = model_combine[model_idx+pianyi].decodeSymbol(&rc);
    }
    else //combine compress
    {
        decoded_matrix_seq[offset] = seq[tmp];
        decoded_matrix_qual[offset] = qual[tmp];
        if(int(decoded_matrix_seq[offset])==0) decoded_matrix_seq[offset]='A';
        else if(int(decoded_matrix_seq[offset])==1) decoded_matrix_seq[offset]='T';
        else if(int(decoded_matrix_seq[offset])==2) decoded_matrix_seq[offset]='C';
        else if(int(decoded_matrix_seq[offset])==3) decoded_matrix_seq[offset]='G';
        else if(int(decoded_matrix_seq[offset])==4) decoded_matrix_seq[offset]='N';
    }
     }
  }
    for(int tmp = 0;tmp<sourceLen;tmp++)
    {
        decoded_matrix_qual[tmp] = decoded_matrix_qual[tmp] + min_1;
    }

    rc.FinishDecode();


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
