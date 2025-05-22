#include "HashRefIndex.h"
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

char *IHashRefIndex::createShm(uint64_t lsize)
{
    string str(m_paramPtr->m_strRefFile);
	str = str.substr(str.rfind('/')+1);

	int shmid; 
	if ((shmid = shm_open(str.c_str(), O_RDWR, 0)) < 0) {
		shmid = shm_open(str.c_str(), O_CREAT|O_RDWR|O_EXCL, 0644);
	}
	else
	{
		struct stat statbuf;
		fstat(shmid, &statbuf);
		if(statbuf.st_size == m_iRefFileSize)
		{
			arc_stdout(level::INFO, "createHashShm %s has existed", str.c_str());
			return nullptr;
		}
		else
		{
			shm_unlink(str.c_str());
			shmid = shm_open(str.c_str(), O_CREAT|O_RDWR|O_EXCL, 0644);
		}
	}
	if (shmid < 0) 
	{
		arc_stdout(level::ERROR, "shm_open fail");
		return nullptr;
	}

	if(ftruncate(shmid, lsize) == -1)
	{
		arc_stdout(level::ERROR, "shm_open fail");
		return nullptr;
	}
	char *shm = (char*)mmap(0, lsize, PROT_READ|PROT_WRITE, MAP_SHARED, shmid, 0);
	if(shm == MAP_FAILED)
	{
		arc_stdout(level::ERROR, "mmap fail");
		return nullptr;
	}

    return shm;
}

HashRefIndex32::HashRefIndex32(uint64_t filesize):
    IHashRefIndex(filesize)
{
}

HashRefIndex32::~HashRefIndex32()
{
    if(m_bNeedFree)
    {
        if(m_pSeqint) free(m_pSeqint);
        if(m_pSeednum) free(m_pSeednum);
        if(m_pSeedind) free(m_pSeedind);
        if(m_pSeedpos) free(m_pSeedpos);
    }
}


void HashRefIndex32::initMemory()
{
	uint64_t maxnseqint = (m_iRefFileSize>>4)+1;
    m_pSeqint = (uint32_t*) calloc(maxnseqint, sizeof(uint32_t));
    m_pSeednum = (uint32_t*) calloc(m_iMaxnseeds, sizeof(uint32_t));
}


bool HashRefIndex32::createRefIndexShm()//共享区
{
    uint64_t total = 4+m_iSeqint+m_iMaxnseeds*2+m_iTotnpos;
	total <<= 2;
    char *shm = createShm(total);
    if(shm == nullptr)
    {
		arc_stdout(level::ERROR, "createShm fail");
        return false;
    }

	IntTo4Ch(m_iSeedlen, shm);
	shm += sizeof(uint32_t);
	IntTo4Ch(m_iChnum, shm);
	shm += sizeof(uint32_t);
	IntTo4Ch(m_iSeqint, shm);
	shm += sizeof(uint32_t);
	IntTo4Ch(m_iTotnpos, shm);
	shm += sizeof(uint32_t);

	memcpy(shm, m_pSeqint, sizeof(uint32_t)*m_iSeqint);
	shm += sizeof(uint32_t)*m_iSeqint;
	memcpy(shm, m_pSeednum, sizeof(uint32_t)*m_iMaxnseeds);
	shm += sizeof(uint32_t)*m_iMaxnseeds;
	memcpy(shm, m_pSeedind, sizeof(uint32_t)*m_iMaxnseeds);
	shm += sizeof(uint32_t)*m_iMaxnseeds;
	memcpy(shm, m_pSeedpos, sizeof(uint32_t)*m_iTotnpos);

	return true;
}


bool HashRefIndex32::loadRefIndexShm()
{
    m_bNeedFree = false;
    string str(m_paramPtr->m_strRefFile);
	str = str.substr(str.rfind('/')+1);

	int shmid; 
	if ((shmid = shm_open(str.c_str(), O_RDONLY, 0)) < 0) return false;
	struct stat statbuf;
	fstat(shmid, &statbuf);
	unsigned char *shm = (unsigned char*)mmap(0, 16, PROT_READ, MAP_SHARED, shmid, 0);
	
	m_iSeedlen = DECODE_INT(shm);
	shm += sizeof(uint32_t);
	m_iChnum = DECODE_INT(shm);
	shm += sizeof(uint32_t);
	m_iSeqint = DECODE_INT(shm);
	shm += sizeof(uint32_t);
	m_iTotnpos = DECODE_INT(shm);
	shm += sizeof(uint32_t);

    m_iMaxnseeds = g_mask[m_iSeedlen]+1;
	uint64_t total = 4+m_iSeqint+m_iMaxnseeds*2+m_iTotnpos;
	total <<= 2;

	if(total != statbuf.st_size)
	{
		arc_stdout(level::ERROR, "/dev/shm/ %s is wrong file, please delete", str.c_str());
		return false;
	}
	shm = (unsigned char*)mmap(0, total, PROT_READ, MAP_SHARED, shmid, 0);
	shm += 16;
	m_pSeqint = (uint32_t*)shm;
	shm += sizeof(uint32_t)*m_iSeqint;
	m_pSeednum = (uint32_t*)shm;
	shm += sizeof(uint32_t)*m_iMaxnseeds;
	m_pSeedind = (uint32_t*)shm;
	shm += sizeof(uint32_t)*m_iMaxnseeds;
	m_pSeedpos = (uint32_t*)shm;
	// cout << "m_pSeedpos:" << *m_pSeedpos << endl;

	return true;
}


bool HashRefIndex32::loadRefIndexFile()
{
    m_bNeedFree = true;
    char refpath[1000];
    sprintf(refpath, "%s.hash", m_paramPtr->m_strRefFile.c_str());

    FILE *fn = fopen(refpath,"r");
    if(fn == NULL){
		arc_stdout(level::ERROR, "read info file failed");
    } else{
        fread(&m_iSeedlen, sizeof(uint32_t), 1, fn);//从文件里面读取
        fread(&m_iChnum, sizeof(uint32_t), 1, fn);
        fread(&m_iSeqint, sizeof(uint32_t), 1, fn);
        fread(&m_iTotnpos, sizeof(uint32_t), 1, fn);
        
        m_iMaxnseeds = g_mask[m_iSeedlen]+1;
        m_pSeqint = (uint32_t*) calloc(m_iSeqint, sizeof(uint32_t)); 
        fread(m_pSeqint, sizeof(uint32_t), m_iSeqint, fn);//读取m_iSeqint个字符

		if(m_paramPtr->m_actionType == ACTIONTYPE_DOENCODE)
		{
			m_pSeednum = (uint32_t*) calloc(m_iMaxnseeds, sizeof(uint32_t)); // the number of occurences of each seed
        	m_pSeedind = (uint32_t*) calloc(m_iMaxnseeds, sizeof(uint32_t)); // the indices of positions of each seed
        	m_pSeedpos = (uint32_t*) calloc(m_iTotnpos, sizeof(uint32_t)); // all positions ordering by seeds

			fread(m_pSeednum, sizeof(uint32_t), m_iMaxnseeds, fn);
			fread(m_pSeedind, sizeof(uint32_t), m_iMaxnseeds, fn);
			fread(m_pSeedpos, sizeof(uint32_t), m_iTotnpos, fn);

			createRefIndexShm();
		}
    }
    
    fclose(fn);
	return true;
}


void HashRefIndex32::setSeqint(uint64_t Seqint, uint8_t chval)
{
	m_pSeqint[Seqint] = (m_pSeqint[Seqint]<<2)|chval;
}


void HashRefIndex32::setEndSeqint(uint64_t Chnum, uint64_t Seqint)
{
	int left = Chnum%m_ihalf;
	if(left != 0)
	{
		m_pSeqint[Seqint] = m_pSeqint[Seqint]<<(2*(m_ihalf-left));
	}
    m_iChnum = Chnum;
	m_iSeqint = Seqint+1;
}


void HashRefIndex32::setSeednum(uint32_t seedseq)
{
    if(m_pSeednum[seedseq] < m_ilimit){
        m_pSeednum[seedseq]++;
    }
}

/**
 * @brief  分配内存，并给m_pSeedind赋值
 * @param  mask: 掩码
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:17:53
 */
void HashRefIndex32::setSeedind(uint64_t mask)
{
	m_pSeedind = (uint32_t*) calloc(m_iMaxnseeds, sizeof(uint32_t));
    for(uint32_t i = 0; i < mask; i++){
		if(m_pSeednum[i] >= m_ilimit)
			m_pSeednum[i] = 0;
        m_pSeedind[i+1] = m_pSeedind[i]+m_pSeednum[i];
        m_iTotnpos += m_pSeednum[i];
    }
	if(m_pSeednum[mask] >= m_ilimit)
		m_pSeednum[mask] = 0;
	m_iTotnpos += m_pSeednum[mask];

    m_pSeedpos = (uint32_t*) calloc(m_iTotnpos, sizeof(uint32_t)); 
    m_pSeedindshift = (uint32_t*) calloc(m_iMaxnseeds, sizeof(uint32_t));
}

/**
 * @brief  给m_pSeedpos数组赋值
 * @param  seedseq: m_pSeednum下标
 * @param  Chnum: 参考序列的碱基个数
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:18:37
 */
void HashRefIndex32::setSeedpos(uint32_t seedseq, uint64_t Chnum)
{
    if(m_pSeednum[seedseq] > 0){
        m_pSeedpos[m_pSeedind[seedseq]+m_pSeedindshift[seedseq]] = Chnum - m_ishort;
        m_pSeedindshift[seedseq]++;
    }
}

/**
 * @brief  写索引文件到磁盘
 * @retval 成功为true,失败false
 * @author  zhaozijian 2019/9/20 上午11:19:25
 */
bool HashRefIndex32::writeIndexFile()
{
    free(m_pSeedindshift);
    char hashpath[PATH_LEN]={0};
    sprintf(hashpath, "%s.hash", m_paramPtr->m_strRefFile.c_str());
    FILE *fhash = fopen(hashpath,"wb+");
    if(fhash == NULL){
		arc_stdout(level::ERROR, "build ref-index file failed");
    }
    else{
        fwrite(&m_iSeedlen, sizeof(uint32_t), 1, fhash);
        fwrite(&m_iChnum, sizeof(uint32_t), 1, fhash);
        fwrite(&m_iSeqint, sizeof(uint32_t), 1, fhash);
        fwrite(&m_iTotnpos, sizeof(uint32_t), 1, fhash);

        fwrite(m_pSeqint, sizeof(uint32_t), m_iSeqint, fhash);
        fwrite(m_pSeednum, sizeof(uint32_t), m_iMaxnseeds, fhash);
        fwrite(m_pSeedind, sizeof(uint32_t), m_iMaxnseeds, fhash);
        fwrite(m_pSeedpos, sizeof(uint32_t), m_iTotnpos, fhash);
		fclose(fhash); 

		createRefIndexShm();
		return true;
    }
    return false;
}

/********************************************************************/

HashRefIndex64::HashRefIndex64(uint64_t filesize):
    IHashRefIndex(filesize)   
{
}

HashRefIndex64::~HashRefIndex64()
{
    if(m_bNeedFree)
    {
        if(m_pSeqint) free(m_pSeqint);
        if(m_pSeednum) free(m_pSeednum);
        if(m_pSeedind) free(m_pSeedind);
        if(m_pSeedpos) free(m_pSeedpos);
    }
}

void HashRefIndex64::initMemory()
{
    uint64_t maxnseqint = (m_iRefFileSize>>4)+1;
    m_pSeqint = (uint32_t*) calloc(maxnseqint, sizeof(uint32_t));
    m_pSeednum = (uint64_t*) calloc(m_iMaxnseeds, sizeof(uint64_t));
}

/**
 * @brief  创建内存映射文件并填充内容
 * @retval 成功为true,失败false
 * @author  zhaozijian 2019/9/20 上午11:20:08
 */
bool HashRefIndex64::createRefIndexShm()
{
	return false;
	uint64_t t_a = m_iSeqint +1;
	t_a <<= 2;
	uint64_t t_b = m_iMaxnseeds*2+m_iTotnpos;
	t_b <<= 3;
    uint64_t total = t_a +t_b;

    char *shm = createShm(total);
    if(shm == nullptr)
    {
		arc_stdout(level::ERROR, "createShm fail");
        return false;
    }

	IntTo4Ch(m_iSeedlen, shm);
	shm += sizeof(uint32_t);
	IntTo8Ch(m_iChnum, shm);
	shm += sizeof(uint64_t);
	IntTo8Ch(m_iSeqint, shm);
	shm += sizeof(uint64_t);
	IntTo8Ch(m_iTotnpos, shm);
	shm += sizeof(uint64_t);

	memcpy(shm, m_pSeqint, sizeof(uint32_t)*m_iSeqint);
	shm += sizeof(uint32_t)*m_iSeqint;
	memcpy(shm, m_pSeednum, sizeof(uint64_t)*m_iMaxnseeds);
	shm += sizeof(uint64_t)*m_iMaxnseeds;
	memcpy(shm, m_pSeedind, sizeof(uint64_t)*m_iMaxnseeds);
	shm += sizeof(uint64_t)*m_iMaxnseeds;
	memcpy(shm, m_pSeedpos, sizeof(uint64_t)*m_iTotnpos);

	return true;
}

uint64_t chtoint64(unsigned char *ch)
{
	return ch[0] + (ch[1]<<8) + (ch[2]<<16) + (ch[3]<<24) +
        ((uint64_t)ch[4]<<32) + ((uint64_t)ch[5]<<40) + 
		((uint64_t)ch[6]<<48) + ((uint64_t)ch[7]<<56);
}

/**
 * @brief  从内存加载索引文件
 * @retval 成功为true,失败false
 * @author  zhaozijian 2019/9/20 上午11:20:47
 */
bool HashRefIndex64::loadRefIndexShm()
{
	return false;
    m_bNeedFree = false;
    string str(m_paramPtr->m_strRefFile);
	str = str.substr(str.rfind('/')+1);

	int shmid; 
	if ((shmid = shm_open(str.c_str(), O_RDONLY, 0)) < 0) return false;
	struct stat statbuf;
	fstat(shmid, &statbuf);
	unsigned char *shm = (unsigned char*)mmap(0, 28, PROT_READ, MAP_SHARED, shmid, 0);
	
	m_iSeedlen = DECODE_INT(shm);
	shm += sizeof(uint32_t);
	m_iChnum = chtoint64(shm);
	shm += sizeof(uint64_t);
	m_iSeqint = chtoint64(shm);
	shm += sizeof(uint64_t);
	m_iTotnpos = chtoint64(shm);
	shm += sizeof(uint64_t);

    m_iMaxnseeds = g_mask[m_iSeedlen]+1;
	uint64_t t_a = m_iSeqint +1;
	t_a <<= 2;
	uint64_t t_b = m_iMaxnseeds*2+m_iTotnpos;
	t_b <<= 3;
    uint64_t total = t_a +t_b;
	if(total != statbuf.st_size)
	{
		arc_stdout(level::ERROR, "/dev/shm/ %s is wrong file, please delete", str.c_str());
		return false;
	}	
	shm = (unsigned char*)mmap(0, total, PROT_READ, MAP_SHARED, shmid, 0);
	shm += 28;
	m_pSeqint = (uint32_t*)shm;
	shm += sizeof(uint32_t)*m_iSeqint;
	m_pSeednum = (uint64_t*)shm;
	shm += sizeof(uint64_t)*m_iMaxnseeds;
	m_pSeedind = (uint64_t*)shm;
	shm += sizeof(uint64_t)*m_iMaxnseeds;
	m_pSeedpos = (uint64_t*)shm;

	return true;
}

/**
 * @brief  从硬盘加载索引文件
 * @retval 成功为true,失败false
 * @author  zhaozijian 2019/9/20 上午11:21:04
 */
bool HashRefIndex64::loadRefIndexFile()
{
    m_bNeedFree = true;
    char refpath[1000];
    sprintf(refpath, "%s.hash", m_paramPtr->m_strRefFile.c_str());

    FILE *fn = fopen(refpath,"r");
    if(fn == NULL){
		arc_stdout(level::ERROR, "read info file failed");
    } else{
        fread(&m_iSeedlen, sizeof(uint32_t), 1, fn);
        fread(&m_iChnum, sizeof(uint64_t), 1, fn);
        fread(&m_iSeqint, sizeof(uint64_t), 1, fn);
        fread(&m_iTotnpos, sizeof(uint64_t), 1, fn);
        
        m_iMaxnseeds = g_mask[m_iSeedlen]+1;
        m_pSeqint = (uint32_t*) calloc(m_iSeqint, sizeof(uint32_t)); 
        fread(m_pSeqint, sizeof(uint32_t), m_iSeqint, fn);

		if(m_paramPtr->m_actionType == ACTIONTYPE_DOENCODE)
		{
			m_pSeednum = (uint64_t*) calloc(m_iMaxnseeds, sizeof(uint64_t)); // the number of occurences of each seed
        	m_pSeedind = (uint64_t*) calloc(m_iMaxnseeds, sizeof(uint64_t)); // the indices of positions of each seed
        	m_pSeedpos = (uint64_t*) calloc(m_iTotnpos, sizeof(uint64_t)); // all positions ordering by seeds

			fread(m_pSeednum, sizeof(uint64_t), m_iMaxnseeds, fn);
			fread(m_pSeedind, sizeof(uint64_t), m_iMaxnseeds, fn);
			fread(m_pSeedpos, sizeof(uint64_t), m_iTotnpos, fn);

			//createRefIndexShm();
		}
    }
    
    fclose(fn);
	return true;
}

/**
 * @brief  给m_pSeqint数组赋值
 * @param  Seqint: m_pSeqint数组下标
 * @param  chval: 碱基对应的bit值
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:15:14
 */
void HashRefIndex64::setSeqint(uint64_t Seqint, uint8_t chval)
{
	m_pSeqint[Seqint] = (m_pSeqint[Seqint]<<2)|chval;
}

/**
 * @brief  给m_pSeqint数组结尾赋值
 * @param  Chnum: 参考序列的碱基个数
 * @param  Seqint: m_pSeqint数组下标
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:16:36
 */
void HashRefIndex64::setEndSeqint(uint64_t Chnum, uint64_t Seqint)
{
    int left = Chnum%m_ihalf;
	if(left != 0)
	{
		m_pSeqint[Seqint] = m_pSeqint[Seqint]<<(2*(m_ihalf-left));
	}
    m_iChnum = Chnum;
	m_iSeqint = Seqint+1;
}

/**
 * @brief  给m_pSeednum数组赋值
 * @param  seedseq: m_pSeednum数组下标
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:17:27
 */
void HashRefIndex64::setSeednum(uint32_t seedseq)
{
    if(m_pSeednum[seedseq] < m_ilimit){
        m_pSeednum[seedseq]++;
    }
}

/**
 * @brief  分配内存，并给m_pSeedind赋值
 * @param  mask: 掩码
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:17:53
 */
void HashRefIndex64::setSeedind(uint64_t mask)
{
	m_pSeedind = (uint64_t*) calloc(m_iMaxnseeds, sizeof(uint64_t));
    for(uint32_t i = 0; i < mask; i++){
		if(m_pSeednum[i] >= m_ilimit)
			m_pSeednum[i] = 0;
        m_pSeedind[i+1] = m_pSeedind[i]+m_pSeednum[i];
        m_iTotnpos += m_pSeednum[i];
    }
	if(m_pSeednum[mask] >= m_ilimit)
		m_pSeednum[mask] = 0;
	m_iTotnpos += m_pSeednum[mask];

    m_pSeedpos = (uint64_t*) calloc(m_iTotnpos, sizeof(uint64_t));
    m_pSeedindshift = (uint64_t*) calloc(m_iMaxnseeds, sizeof(uint64_t));
}

/**
 * @brief  给m_pSeedpos数组赋值
 * @param  seedseq: m_pSeednum下标
 * @param  Chnum: 参考序列的碱基个数
 * @retval None
 * @author  zhaozijian 2019/9/20 上午11:18:37
 */
void HashRefIndex64::setSeedpos(uint32_t seedseq, uint64_t Chnum)
{
    if(m_pSeednum[seedseq] > 0){
        m_pSeedpos[m_pSeedind[seedseq]+m_pSeedindshift[seedseq]] = Chnum - m_ishort;
        m_pSeedindshift[seedseq]++;
    }
}

/**
 * @brief  写索引文件到磁盘
 * @retval 成功为true,失败false
 * @author  zhaozijian 2019/9/20 上午11:19:25
 */
bool HashRefIndex64::writeIndexFile()
{
    free(m_pSeedindshift);
    char hashpath[PATH_LEN]={0};
    sprintf(hashpath, "%s.hash", m_paramPtr->m_strRefFile.c_str());
    FILE *fhash = fopen(hashpath,"wb+");
    if(fhash == NULL){
		arc_stdout(level::ERROR, "build ref-index file failed ");
    }
    else{
        fwrite(&m_iSeedlen, sizeof(uint32_t), 1, fhash);
        fwrite(&m_iChnum, sizeof(uint64_t), 1, fhash);
        fwrite(&m_iSeqint, sizeof(uint64_t), 1, fhash);
        fwrite(&m_iTotnpos, sizeof(uint64_t), 1, fhash);

        fwrite(m_pSeqint, sizeof(uint32_t), m_iSeqint, fhash);
        fwrite(m_pSeednum, sizeof(uint64_t), m_iMaxnseeds, fhash);
        fwrite(m_pSeedind, sizeof(uint64_t), m_iMaxnseeds, fhash);
        fwrite(m_pSeedpos, sizeof(uint64_t), m_iTotnpos, fhash);
		fclose(fhash); 

		createRefIndexShm();
		return true;
    }
    return false;
}