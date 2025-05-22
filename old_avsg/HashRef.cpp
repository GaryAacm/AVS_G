#include "HashRef.h"
#include "AlignInfo.h"
#include <math.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "Md5.h"

uint16_t g_mismatch_count[65536];

HashRef::HashRef()
{
	uint64_t filesize = getFileSize(m_paramPtr->m_strRefFile.c_str());

	int l = filesize>>30;
	if(l > 4)
	{
		m_hashIndexPtr = new HashRefIndex64(filesize);
	}
	else
	{
		m_hashIndexPtr = new HashRefIndex32(filesize);
	}
}

HashRef::~HashRef()
{
	freeRefIndex();
}

/**
 * @brief  判断参考序列的索引文件是否存在
 * @retval 存在返回true,否则false
 * @author  zhaozijian 2019/9/4 下午2:25:41
 */
bool HashRef::isRefIndexExist()
{
	char hashpath[PATH_LEN]={0};
    sprintf(hashpath, "%s.hash", m_paramPtr->m_strRefFile.c_str());

    char path[PATH_LEN] = {0};
    sprintf(path, "%s.md5", m_paramPtr->m_strRefFile.c_str());

	if(isFileExist(hashpath) && isFileExist(path))
	{
		return true;
	}
	return false;
}

/**
 * @brief  构建参考序列的hash索引文件
 * @retval 构建成功返回true,否则false
 * @author  zhaozijian 2019/9/4 下午2:28:19
 */
bool HashRef::buildRefIndex()
{
	if(isRefIndexExist())
	{
		return true;
	}

	FILE *f = fopen(m_paramPtr->m_strRefFile.c_str(), "r");
    if(f == NULL){
		arc_stdout(level::ERROR, "read reference file failed");
    }

	uint32_t iseed = 0;
    uint32_t nchr = 0;
    uint32_t seedseq = 0;
	uint64_t mask = g_mask[SEEDLEN];
	uint64_t Seqint = 0;
	uint64_t Chnum = 0;
	int half = m_hashIndexPtr->getHalf();
	m_hashIndexPtr->initMemory();

    char* seq;
    size_t seqlen = 0;
    ssize_t linelength;
    uint8_t chval = 4;
    while((linelength = getline(&seq, &seqlen, f)) != -1){
        if(*seq == '>'){
            nchr++;
        } else if(nchr > 0){
            for(int i = 0; i < strlen(seq)-1; i++){
                Seqint = Chnum/half;
                Chnum++;
				chval = hash_std_table[seq[i]] & 3;
                seedseq = mask & ((seedseq << 2)|chval);

				m_hashIndexPtr->setSeqint(Seqint, chval);

                if((seq[i] == 'N')||(seq[i] == 'n')){
                        iseed = 0;
                } else {
                    iseed ++;
                    if(iseed == SEEDLEN){
                        if(Chnum % INTERVAL == 0){
							m_hashIndexPtr->setSeednum(seedseq);
                        }
                        iseed--;
                    }
                }
                	
            }
        }
    }
    fclose(f);
    
    if(nchr == 0){
		arc_stdout(level::ERROR, "read reference file failed");
    }

	m_hashIndexPtr->setEndSeqint(Chnum, Seqint);
	m_hashIndexPtr->setSeedind(mask);
    

    f = fopen(m_paramPtr->m_strRefFile.c_str(),"r");
    Chnum = 0;
    seedseq = 0;
    iseed = 0;
    int i_tmp = SEEDLEN - 1;
    while((linelength = getline(&seq, &seqlen, f)) != -1){
        if(*seq != '>'){
            for(int i = 0; i < strlen(seq)-1; i++){                   
                Chnum++;
				chval = hash_std_table[seq[i]] & 3;
                seedseq = mask & ((seedseq << 2)|chval);
    
                if((seq[i] == 'N')||(seq[i] == 'n')){
                    iseed = 0;
                } else {
                    iseed++;
                    if(iseed == SEEDLEN){
                        if(Chnum % INTERVAL == 0){
							m_hashIndexPtr->setSeedpos(seedseq, Chnum);
                        }
                        iseed--;
                    }
                } 
            }
        }
    }
    fclose(f);

	bool ret = m_hashIndexPtr->writeIndexFile();
	unsigned char md5[16];
    MD5File(m_paramPtr->m_strRefFile.c_str(), md5);
    char path[PATH_LEN] = {0};
    sprintf(path, "%s.md5", m_paramPtr->m_strRefFile.c_str());
    fstream out_s;
    out_s.open(path, std::ios::binary | std::ios::out);
    out_s.write((char*)md5, 16);
    out_s.close();
	return ret;
}

/**
 * @brief  加载索引文件
 * @retval true
 * @author  zhaozijian 2019/9/4 下午2:29:12
 */
bool HashRef::loadRefIndex()
{   
    int k,l;
    for( k = 1 ; k < 65536 ; k++){
        uint16_t mis = 0,tmp = k;
        for( l = 0 ; l < 8 ; l++){
            if((tmp & 0x3) != 0) {
                mis++;
            }
            tmp >>= 2;
        }
        g_mismatch_count[k] = mis;
    }

	bool ret = m_hashIndexPtr->loadRefIndexShm();
	if(!ret)
	{
		m_hashIndexPtr->loadRefIndexFile();
	}
	m_kmerlen = m_hashIndexPtr->getSeedlen();
	m_reflen = m_hashIndexPtr->getChnum();

	m_referencelen = m_hashIndexPtr->getChnum();
	//m_offsetBit = getbitnum(m_referencelen) - 2;
    m_offsetSize =  pow(2, m_offsetBit) - 1;	
	return true;	
}

bool HashRef::freeRefIndex()
{
	RELEASEPTR(m_hashIndexPtr);
	return true;
}

uint32_t HashRef::doSEAlign(AlignInfo &info)
{
	return getHashAlignInfo(info);
}

uint32_t HashRef::doPEAlign(AlignInfo &info1, AlignInfo &info2)
{
	uint32_t mis1 = getHashAlignInfo(info1);
	uint32_t mis2 = getHashAlignInfo(info2);
	return 0;
}

uint32_t HashRef::getHashAlignInfo(AlignInfo &info)
{
	info.m_bRev = false;
	int degenum = getHashSeeds(info);
	if(degenum > info.m_limit) //中间的简并个数超出maxmis，直接返回比对失败
	{
		return info.m_misnum;
	}
	bool match = hashAligner(info);
	
	if(!match && info.m_seedextendcount < MAX_ALIGN_POSITION)//反向互补碱基
	{
		info.m_bRev = true;
		for(int j=info.m_seqlen-1; j >=0; j--){
			if(info.m_seqarry[j] > 3)
			{
				info.m_revseqarry.push_back(info.m_seqarry[j]);
			}
			else
			{
				info.m_revseqarry.push_back(info.m_seqarry[j]^0x3);
			}
		}

		getHashSeeds(info);
		match = hashAligner(info);

		if(m_bslice)
		{
			if(!match && info.m_seedextendcount < MAX_ALIGN_POSITION)
			{
				info.m_bRev = false;
				match = hashAlignerShortPart(info);
			}
			if(!match && info.m_seedextendcount < MAX_ALIGN_POSITION)
			{
				info.m_bRev = true;
				match = hashAlignerShortPart(info);
			}
		}
	}

	return info.m_misnum;
}

int HashRef::getHashSeeds(AlignInfo &info)//找到简并碱基
{
	Memory<uint8_t> &seqarry = info.m_bRev ? info.m_revseqarry : info.m_seqarry;
	Memory<uint32_t> &seqbase = info.m_bRev ? info.m_revseqbase: info.m_seqbase;
	Memory<uint32_t> &fqseed = info.m_bRev ? info.m_revfqseed: info.m_fqseed;
	//uint8_t chval = 5;
	int degenum = 0, ifq = 0;
	uint32_t kmer = 0;
	for(int i = 0; i < info.m_seqlen; i++){
		if(seqarry[i] > 3)
		{
			degenum++;
		}
		//chval = seqarry[i] & 3;															
		ifq = i>>4;	
		//uint32_t wj2 = seqbase[ifq];//wj
		seqbase[ifq] = (seqbase[ifq]<<2) | seqarry[i]; 													
		kmer = g_mask[m_kmerlen] & ((kmer << 2)|seqarry[i]);
		//int wj3 = 0;
		if(i >= (m_kmerlen-1)){
			fqseed.push_back(kmer);//wj 求k-mer
		}
	}
	
	if(info.m_seqlen&0xf)
	{
		int len = (info.m_seqlen - 1) >> 4;
		seqbase[len] <<= (32 - ((info.m_seqlen&0xf) << 1));
	}

	return degenum;
}

//寻找seed个数最小的pos
bool HashRef::findHashSeeds(Memory<uint32_t> &fqseed, int start, int end, int maxdepth, hashseed &hseed, bool bfirst)
{
	int i,minpos = SEED_MAX_DEPTH;
	uint64_t cnt = 0;
	for(i = start; i <= end; i+=2)
	{
		cnt = m_hashIndexPtr->getSeednum(fqseed[i]);//找到全部匹配上的数据
		if(cnt && cnt < maxdepth){
			if (cnt < minpos)
			{
				minpos = cnt;
				hseed.seqpos = i;
				hseed.fqseed = fqseed[i];
				hseed.seednum = cnt;
				if(bfirst) break;
			}
		}
	}
	if(minpos == SEED_MAX_DEPTH)
		return false;
	return true;
}

bool HashRef::gaplessHashAlignPositions(uint64_t refPos, AlignInfo &info)
{
	int i, j, nmis = 0, seqlen=info.m_seqlen;
	uint32_t mode = (refPos -1) & 0xf, shift = mode << 1 , r_shift = 32 - shift,len = (seqlen-1) >> 4;
	uint64_t idx = ((refPos-1) >> 4) , nmask = ~g_mask[mode], seqdot = 0;
	info.m_seedextendcount++;
	Memory<uint8_t> &seqarry = info.m_bRev ? info.m_revseqarry : info.m_seqarry;
	Memory<uint32_t> &seqbaseArry = info.m_bRev ? info.m_revseqbase : info.m_seqbase;
	for(j = 0 ; ( j <= len && nmis <= info.m_limit ) ; j++){
		seqdot = 0;
		uint32_t seqint_j = m_hashIndexPtr->getSeqint(idx+j);
		if(mode)
		{
			uint32_t seqint_j_1 = m_hashIndexPtr->getSeqint(idx+j+1);
			seqdot = seqbaseArry[j]^ (((seqint_j << shift) & nmask )| ((seqint_j_1 >> r_shift) & g_mask[mode]));
		} else {
			seqdot = seqbaseArry[j]^(seqint_j);
		}
		
		if(j >= (seqlen>>4))
			seqdot &= ~g_mask[16-(seqlen&0xf)];
		
		if(seqdot){
			for( i = 0 ; i < 2 ; i++){
				if((seqdot & 0xffff) != 0) {
					nmis += g_mismatch_count[seqdot & 0xffff];
				}
				seqdot >>= 16;
			}
		}
	}

	if(nmis<= info.m_limit)
	{
		nmis = 0;
		uint8_t maxmis = info.m_limit+1;
		uint32_t start, tmp, refBaseInt, readBaseInt;
		info.m_cigarL.clear();
		info.m_cigarV.clear();
		for(j = 0; j < seqlen; j++)
		{
			start = refPos + j - 1; // refPos start from 1 not 0
			mode = start & 0xf;
			idx = start>>4;
			tmp = m_hashIndexPtr->getSeqint(idx);
			refBaseInt = tmp >> (30 - (mode << 1)) & 3;
			readBaseInt = seqarry[j];

			if(refBaseInt != readBaseInt){
				if(nmis < maxmis)
				{
					nmis++;
				}
				else
				{
					goto END;
				}
				
				info.m_cigarL.push_back(j);
				if(readBaseInt > 3){
					info.m_cigarV.push_back(3);
				} else {
					info.m_cigarV.push_back(basesHash[refBaseInt][readBaseInt]);
				}
			}
		}
		
		if(nmis < maxmis){
			info.m_alignpos = refPos;
			info.m_misnum = nmis;
			return true;
		}
	}

END:
	return false;
}

bool HashRef::gaplessHashAlign(int idx, AlignInfo &info)
{
	int MinMis = info.m_limit + 1;
	uint32_t MinPos = 0;

	bool match = false;
	uint64_t refpos;
	hashseed &hseed = info.m_hashseed[idx];
	for(uint64_t j = 0; j < hseed.seednum && info.m_seedextendcount <= MAX_ALIGN_POSITION; j++)
	{
		refpos = m_hashIndexPtr->getSeedpos(hseed.fqseed, j);
		refpos -= hseed.seqpos;
		if(refpos > 0 && refpos < (m_reflen - info.m_seqlen))
		{					
			match = gaplessHashAlignPositions(refpos, info);
			if(match)
			{
				if (info.m_misnum < MinMis)
				{
					MinMis = info.m_misnum;
					MinPos = info.m_alignpos;
					if (MinMis == 0)
					{
						break;
					}
				}
				// return match;
			}
		}
	}
	if (MinMis < info.m_limit + 1)
	{
		return gaplessHashAlignPositions(MinPos, info);
	}

	return false;
}
// bool HashRef::gaplessHashAlign(int idx, AlignInfo &info)
// {
// 	bool match = false;
// 	uint64_t refpos;
// 	hashseed &hseed = info.m_hashseed[idx];
// 	for(uint64_t j = 0; j < hseed.seednum && info.m_seedextendcount <= MAX_ALIGN_POSITION; j++)
// 	{
// 		refpos = m_hashIndexPtr->getSeedpos(hseed.fqseed, j);
// 		refpos -= hseed.seqpos;
// 		if(refpos > 0 && refpos < (m_reflen - info.m_seqlen))
// 		{					
// 			match = gaplessHashAlignPositions(refpos, info);
// 			if(match)
// 			{
// 				return match;
// 			}
// 		}
// 	}
// 	return false;
// }


bool HashRef::hashAligner(AlignInfo &info)
{
	Memory<uint32_t> &fqseedarry = info.m_bRev ? info.m_revfqseed : info.m_fqseed;
	for(int i = 0 ; i < 2 ; i++)
	{
		hashseed &hseed = info.m_hashseed[i];
		bool find = findHashSeeds(fqseedarry, i, (info.m_seqlen-m_kmerlen), SEED_MAX_DEPTH, hseed, false);
		if(find)
		{
			bool ret = gaplessHashAlign(i, info);
			if(ret)
			{
				return ret;
			}
		}
	}
	
	return false;
}

bool HashRef::hashAlignerShortPart(AlignInfo &info)
{
	int seqlen = info.m_seqlen;
	int iter = (seqlen > 75) ? 4 : ((seqlen < 45) ? 2 : 3);
	int isbg = (seqlen > iter*m_kmerlen) ? (seqlen/iter-m_kmerlen) : 0;
	Memory<uint32_t> &fqseedarry = info.m_bRev ? info.m_revfqseed : info.m_fqseed;
	bool find = false;
	for(int i = 0 ; i < iter ; i++){
		for(int j = 2 ; j < 4 ; j++){
			hashseed &hseed = info.m_hashseed[j];
			find = findHashSeeds(fqseedarry, i*seqlen/iter+j-2, (i+1)*seqlen/iter - m_kmerlen, PART_SEED_MAX_DEPTH,hseed, false);
			if(find && hseed.seqpos != info.m_hashseed[j-2].seqpos)
			{
				bool ret = gaplessHashAlign(j, info);
				if(ret)
				{
					return ret;
				}
			}
		}
	}

	if(!find){
		for(int i = 0 ; i < iter ; i++){
			for(int j = 2 ; j < 4 ; j++){
				int end = (((i+1)*seqlen/iter-8) > (seqlen - m_kmerlen) ) ? (seqlen - m_kmerlen) : ((i+1)*seqlen/iter-8);
				hashseed &hseed = info.m_hashseed[j];
				find = findHashSeeds(fqseedarry, isbg+i*seqlen/iter+j-2, end, PART_SEED_MAX_DEPTH, hseed, true);
				if(find && hseed.seqpos != info.m_hashseed[j-2].seqpos)
				{
					bool ret = gaplessHashAlign(j, info);
					if(ret)
					{
						return ret;
					}
				}
			}
		}
	}
	return false;
}

int HashRef::query(uint64_t alignPos, int seqlen, uint8_t *refseq)
{
    for(int i = 0; i < seqlen; i++){
        uint64_t start = alignPos + i - 1; 
        uint32_t mode = start & 0xf;
        uint32_t idx = start>>4;
		uint32_t tmp = m_hashIndexPtr->getSeqint(idx);
        refseq[i] = tmp >> (30 - (mode << 1)) & 3;
    }
	return 0;
}