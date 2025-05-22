#include "Context.h"
#include "HashRef.h"
#include "ArcFile.h"
#include "FastqRead.h"
#include "ThreadPool.h"
#include "IWorker.h"
#include "DecodeWorker.h"
#include "DecodeSEWorker.h"
#include "DecodePipeWorker.h"
#include "DecodePEWorker.h"
#include "ExtractWorker.h"
#include "log/progress_bar.h"
#include <fcntl.h> 

Context::Context()
{
    m_paramPtr = Param::GetInstance();
}

Context::~Context()
{
    RELEASEPTR(m_refPtr);
}

void Context::initEncodeContext()
{
    createRefPtr();
    doLoadIndex();
    doReadAndEncode();
}

void Context::initDecodeContext()
{
    createRefPtr();
    doLoadIndex();
    doDecode();
}

void Context::initExtractContext()
{
    if(m_paramPtr->m_extracttype == EXTRACT_LIST)
    {
        printf("NO.\t--------\tFile\n");
        for(int i=0;i<m_paramPtr->m_vecFqFile.size();i++)
        {
            printf("%5d\t--------\t%s\n", i, m_paramPtr->m_vecFqFile[i].c_str());
        }
        return;
    }
    createRefPtr();
    doLoadIndex();
    doExtract();
}

void Context::createRefPtr()
{
    if(m_paramPtr->m_alignType == ALIGNTYPE_HASH)
    {
        m_refPtr = new HashRef();
    }
}

bool Context::doBuildIndex()
{
    createRefPtr();
    if(m_refPtr)
    {
        return m_refPtr->buildRefIndex();
    }
    return false;
}

bool Context::doLoadIndex()
{
    if(m_refPtr)
    {
        m_refPtr->buildRefIndex();
        return m_refPtr->loadRefIndex();
    }
    return false;
}

void Context::doReadAndEncode()
{
    ProgressBar probar;
    int totalcount = m_paramPtr->m_totalfilesize/m_paramPtr->m_iBlockSize + 1;
    probar.setTotal(totalcount);

    ArcFile arcfile;
    FastqRead *readPtr = new FastqRead(m_refPtr);
    readPtr->setProcessBar(&probar); //读取也计算进度
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        readPtr->ReadPEFastq();
    }
    else
    {
        readPtr->ReadFastq();
    }
    
    int current = readPtr->getBlockNum()/3;
    while (true)
    {
        if(readPtr->isFinish())
        {
            break;
        }

        if(current < readPtr->getProcessCount())
        {
            probar.update(readPtr->getProcessCount());
        }
        sleep(1);
    }

    RELEASEPTR(readPtr);
    uint64_t arcsz = arcfile.writeFileInfo();
    probar.end();

    Profile::getInstance().outputCompressResult(arcsz);
    Profile::getInstance().outputElapsed();
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        Profile::getInstance().outputPEAlignmentMisResult(m_paramPtr->m_iBlocknum);
    }
    else
    {
        Profile::getInstance().outputSEAlignmentMisResult(m_paramPtr->m_iBlocknum);
    }
}

void Context::doDecode()
{
    if(m_paramPtr->m_PipeOut)
    {
        pipeOutDecode();
        return;
    }
	
    ThreadPool *m_thPool = new ThreadPool();
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        for(uint32_t i=0;i<m_paramPtr->m_iThreadnum;i++)
        {
            IWorker *ptr = new DecodePEWorker(i, m_refPtr);
            ptr->init();
            m_thPool->addJob(ptr);
        }
    }
    else
    {
        for(uint32_t i=0;i<m_paramPtr->m_iThreadnum;i++)
        {
            IWorker *ptr = new DecodeSEWorker(i, m_refPtr);
            ptr->init();
            m_thPool->addJob(ptr);
        }
    }

    ProgressBar probar;
    probar.setTotal(m_paramPtr->m_iBlocknum+5);
    while (true)
    {
        if(m_thPool->isAllFinish())
        {
            break;
        }
        probar.update(m_thPool->getProcessCount());
        sleep(1);
    }
    probar.end();
    RELEASEPTR(m_thPool);
    Profile::getInstance().outputElapsed();
}

void Context::pipeOutDecode()
{
    ThreadPool *m_thPool = new ThreadPool();

    for(uint32_t i=0;i<m_paramPtr->m_iThreadnum;i++)
    {
        IWorker *ptr = new DecodePipeWorker(i, m_refPtr);
        ptr->init();
        m_thPool->addJob(ptr);
    }

    int blockidx = 0;
    int oldfileidx = -1;
    uint64_t offset[2] ={0,0};
    while (1)
    {
        if(m_thPool->isAllFinish() && m_paramPtr->m_memBufPoolPtr->getFullbufNum() == 0)
        {
            break;
        }

        MemBuf *ptr = m_paramPtr->m_memBufPoolPtr->getIdxFullbuf(blockidx);
        if(ptr == nullptr) continue;
        blockidx++;
        fwrite(ptr->m_bufptr[0], ptr->m_len[0], 1, stdout);
        m_paramPtr->m_memBufPoolPtr->addEmptybuf(ptr);
    }

    RELEASEPTR(m_thPool);
}

void Context::doExtract()
{
    ThreadPool *m_thPool = new ThreadPool();

    for(uint32_t i=0;i<m_paramPtr->m_iThreadnum;i++)
    {
        IWorker *ptr = new ExtractWorker(i, m_refPtr);
        ptr->init();
        m_thPool->addJob(ptr);
    }

    ProgressBar probar;
    probar.setTotal(m_paramPtr->m_iBlocknum+5);

    int blockidx = 0;
    int oldfileidx = -1;
    uint64_t offset[2] ={0,0};
    while (1)
    {
        probar.update(m_thPool->getProcessCount());
        if(m_thPool->isAllFinish() && m_paramPtr->m_memBufPoolPtr->getFullbufNum() == 0)
        {
            break;
        }

        MemBuf *ptr = m_paramPtr->m_memBufPoolPtr->getIdxFullbuf(blockidx);
        if(ptr == nullptr) continue;
        blockidx++;
        switch (m_paramPtr->m_arcType)
        {
            case ARCTYPE_PE:
                pwrite(m_paramPtr->m_outfd[0], ptr->m_bufptr[0], ptr->m_len[0], offset[0]);
                pwrite(m_paramPtr->m_outfd[1], ptr->m_bufptr[1], ptr->m_len[1], offset[1]);
                offset[0] += ptr->m_len[0];
                offset[1] += ptr->m_len[1];
                break;
            case ARCTYPE_SE:
                pwrite(m_paramPtr->m_outfd[0], ptr->m_bufptr[0], ptr->m_len[0], offset[0]);
                offset[0] += ptr->m_len[0];
                break;
            case ARCTYPE_MULTI:
                if(oldfileidx != ptr->m_fileidx)
                {
                    close(m_paramPtr->m_outfd[0]);
                    m_paramPtr->m_outfd[0] = open(m_paramPtr->m_vecFqFile[ptr->m_fileidx].c_str(), O_RDWR|O_CREAT, 0664);
                    oldfileidx = ptr->m_fileidx;
                    offset[0] = 0;
                }
                pwrite(m_paramPtr->m_outfd[0], ptr->m_bufptr[0], ptr->m_len[0], offset[0]);
                offset[0] += ptr->m_len[0];
                break;
            default:
                break;
        }

        m_paramPtr->m_memBufPoolPtr->addEmptybuf(ptr);
    }

    probar.end();
    RELEASEPTR(m_thPool);
    Profile::getInstance().outputElapsed();
}

