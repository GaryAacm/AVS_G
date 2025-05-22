#include "DecodeWorker.h"

DecodeWorker::DecodeWorker(int num, IRef* refptr):
    IWorker(num, refptr)
{
    m_file.open(m_paramPtr->m_strArcFile, std::ios::binary | std::ios::in);
    m_inbuf.mem_malloc(m_paramPtr->m_iBlockSize);
    m_count = m_paramPtr->m_vecBlockMeta.size();
}

DecodeWorker::~DecodeWorker()
{
}

void DecodeWorker::doJob()
{
    auto timer1 = m_profile->getTimer();

    for(uint32_t i=0;i<m_count;i++)
    {
        if(i % m_paramPtr->m_iThreadnum != m_num)
        {
            continue;
        }

        timer1->start("ReadTime");
        elementClear();

        readData(i);
        timer1->stop();

        timer1->start("DecodeTime");
        decodeData();
        timer1->stop();

        timer1->start("RecoverTime");
        uint32_t len = 0;
        if(m_blockParamElPtr->m_onech)
        {
            len = recoverDataNoComment();
        }
        else
        {
            len = recoverDataWithComment();
        }
        timer1->stop();

        timer1->start("WriteTime");
        outPutData(len);
        timer1->stop();
        m_processcount++;
    }

    m_file.close();
}

void DecodeWorker::readData(uint32_t idx)
{
    BlockMetaInfo &info = m_paramPtr->m_vecBlockMeta[idx];
    m_file.seekg(info.compressOffset);
    m_file.read(m_inbuf.getdata(), info.compressSize);

    m_readfileidx = info.fileidx;
    m_blocksize[0] = info.originalSize1;
    m_blocksize[1] = info.originalSize2;
    m_blockoffset[0] = info.originalOffset1;
    m_blockoffset[1] = info.originalOffset2;

}

void DecodeWorker::pe_reorder_process() 
{
    {
        char * tmp = new char[m_seqElPtr->m_rawBuf.size()]; // tmp can be used by seq and qual
        std::vector<char *> tmp_start_ptr(m_seqLenElPtr->m_LenArry.size() + 1);
        tmp_start_ptr[0] = tmp;
        for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); ++i)
            tmp_start_ptr[i + 1] = tmp_start_ptr[i] + m_seqLenElPtr->m_LenArry[i];

        {
            std::memcpy(tmp, m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size());
            char *dest = m_seqElPtr->m_rawBuf.getdata();
            size_t file_1_pos = 0, file_2_pos = 0;
            while((file_1_pos < m_seqLenElPtr->m_LenArry.size()) &&
                  (file_2_pos < m_seqLenElPtr->m_LenArry.size())) {
                while (m_refFreeSeqProcessPtr->pe_id[file_1_pos] != 0) file_1_pos++;
                while (m_refFreeSeqProcessPtr->pe_id[file_2_pos] != 1) file_2_pos++;
                std::memcpy(dest, tmp_start_ptr[file_1_pos], m_seqLenElPtr->m_LenArry[file_1_pos]);
                dest += m_seqLenElPtr->m_LenArry[file_1_pos];
                std::memcpy(dest, tmp_start_ptr[file_2_pos], m_seqLenElPtr->m_LenArry[file_2_pos]);
                dest += m_seqLenElPtr->m_LenArry[file_2_pos];
                file_1_pos++;
                file_2_pos++;
            }
        }

        {
            std::memcpy(tmp, m_qualElPtr->m_rawBuf.getdata(), m_qualElPtr->m_rawBuf.size());
            char * dest = m_qualElPtr->m_rawBuf.getdata();

            size_t file_1_pos = 0, file_2_pos = 0;//读取两个文件测序
            while((file_1_pos < m_seqLenElPtr->m_LenArry.size()) &&
                  (file_2_pos < m_seqLenElPtr->m_LenArry.size())) 
            {
                while (m_refFreeSeqProcessPtr->pe_id[file_1_pos] != 0) file_1_pos++;
                while (m_refFreeSeqProcessPtr->pe_id[file_2_pos] != 1) file_2_pos++;

                std::memcpy(dest, tmp_start_ptr[file_1_pos], m_seqLenElPtr->m_LenArry[file_1_pos]);
                dest += m_seqLenElPtr->m_LenArry[file_1_pos];
                std::memcpy(dest, tmp_start_ptr[file_2_pos], m_seqLenElPtr->m_LenArry[file_2_pos]);
                dest += m_seqLenElPtr->m_LenArry[file_2_pos];

                file_1_pos++;
                file_2_pos++;
            }
        }

        delete [] tmp;
    }

    {
        char * tmp = new char[m_nameElPtr->m_rawBuf.size()];
        std::vector<char *> tmp_start_ptr(m_namelenElPtr->m_LenArry.size() + 1);
        tmp_start_ptr[0] = tmp;
        for (size_t i = 0; i < m_namelenElPtr->m_LenArry.size(); ++i)
            tmp_start_ptr[i + 1] = tmp_start_ptr[i] + m_namelenElPtr->m_LenArry[i];

        std::memcpy(tmp, m_nameElPtr->m_rawBuf.getdata(), m_nameElPtr->m_rawBuf.size());
        char * dest = m_nameElPtr->m_rawBuf.getdata();
        size_t file_1_pos = 0, file_2_pos = 0;
        while((file_1_pos < m_namelenElPtr->m_LenArry.size()) &&
              (file_2_pos < m_namelenElPtr->m_LenArry.size())) {
            while (m_refFreeSeqProcessPtr->pe_id[file_1_pos] != 0) file_1_pos++;
            while (m_refFreeSeqProcessPtr->pe_id[file_2_pos] != 1) file_2_pos++;
            std::memcpy(dest, tmp_start_ptr[file_1_pos], m_namelenElPtr->m_LenArry[file_1_pos]);
            dest += m_namelenElPtr->m_LenArry[file_1_pos];
            std::memcpy(dest, tmp_start_ptr[file_2_pos], m_namelenElPtr->m_LenArry[file_2_pos]);
            dest += m_namelenElPtr->m_LenArry[file_2_pos];
            file_1_pos++;
            file_2_pos++;
        }
        delete [] tmp;
    }

    {
        auto * tmp = new uint32_t [m_seqLenElPtr->m_LenArry.size()];
        size_t file_1_pos = 0, file_2_pos = 0, i = 0;
        std::memcpy(tmp, m_seqLenElPtr->m_LenArry.getdata(), sizeof(uint32_t) * m_seqLenElPtr->m_LenArry.size());
        while((file_1_pos < m_seqLenElPtr->m_LenArry.size()) &&
              (file_2_pos < m_seqLenElPtr->m_LenArry.size())) {
            while (m_refFreeSeqProcessPtr->pe_id[file_1_pos] != 0) file_1_pos++;
            while (m_refFreeSeqProcessPtr->pe_id[file_2_pos] != 1) file_2_pos++;
            m_seqLenElPtr->m_LenArry[i++] = tmp[file_1_pos++];
            m_seqLenElPtr->m_LenArry[i++] = tmp[file_2_pos++];
        }

        file_1_pos = file_2_pos = i = 0;
        std::memcpy(tmp, m_namelenElPtr->m_LenArry.getdata(), sizeof(uint32_t) * m_namelenElPtr->m_LenArry.size());
        while((file_1_pos < m_namelenElPtr->m_LenArry.size()) &&
              (file_2_pos < m_namelenElPtr->m_LenArry.size())) {
            while (m_refFreeSeqProcessPtr->pe_id[file_1_pos] != 0) file_1_pos++;
            while (m_refFreeSeqProcessPtr->pe_id[file_2_pos] != 1) file_2_pos++;
            m_namelenElPtr->m_LenArry[i++] = tmp[file_1_pos++];
            m_namelenElPtr->m_LenArry[i++] = tmp[file_2_pos++];
        }
        delete [] tmp;
    }

}

void DecodeWorker::decodeData()
{
    EncapInfo info;
    char *dataPtr = m_inbuf.getdata();
    Encap::getValue(dataPtr, info);
    if(info.val != BLOCK_ENCAP)
    {
        arc_stdout(level::ERROR, "Wrong block id");
    }
    dataPtr += info.len;
    Encap::getValue(dataPtr, info);
    dataPtr += info.len;
    uint32_t out_len = info.val;
    uint32_t readlen = 0;
    FILE *ss = fopen("decompress.txt","a");
    fprintf(ss,"out_len: %u\n",out_len);
    fclose(ss);
    while (readlen < out_len)
    {
        Encap::getValue(dataPtr, info);
        int id = info.val;
        dataPtr += info.len;
        readlen += info.len;
        Encap::getValue(dataPtr, info);
        
        dataPtr += info.len;
        readlen += info.len;
        

        switch (id)
        {
            case BLOCK_PARAM:
                m_dparam.size = info.val;
                m_dparam.pdata = dataPtr;
                break;
            case BLOCK_SEQLEN:
                m_dseqlen.size = info.val;
                m_dseqlen.pdata = dataPtr;
                break;
            case BLOCK_NAMELEN:
                m_dnamelen.size = info.val;
                m_dnamelen.pdata = dataPtr;
                break;
            case BLOCK_NAME:
                m_dname.size = info.val;
                m_dname.pdata = dataPtr;
                break;
            case BLOCK_SEQ:
                m_dseq.size = info.val;
                m_dseq.pdata = dataPtr;
                break;
            case BLOCK_QUAL:
                m_dqual.size = info.val;
                m_dqual.pdata = dataPtr;
                break;
            case BLOCK_ORDER:
                break;
            case BLOCK_COMBINE:
                m_combine.size = info.val;
                m_combine.pdata = dataPtr;
                ss = fopen("decompress.txt","a");
                fprintf(ss,"size :%u\n",m_dseq.size);
                fclose(ss);
            default:
                break;
        }
        dataPtr += info.val;
        readlen += info.val;
    }

    /*FILE *f = fopen("test_seq.txt", "w");
    char *ptr = m_seqElPtr->m_rawBuf.getdata();
    for(int i=0;i<m_seqLenElPtr->m_LenArry.size();i++)
    {
        fwrite(ptr, 1, m_seqLenElPtr->m_LenArry[i], f);
        fwrite("\n", 1, 1, f);
        ptr += m_seqLenElPtr->m_LenArry[i];
    }
    fclose(f);*/

    auto timer1 = m_profile->getTimer();
    m_paramProcessPtr->decompress(m_dparam.size, m_dparam.pdata);

    timer1->start("SeqLenTime");
    m_seqlenProcessPtr->decompress(m_dseqlen.size, m_dseqlen.pdata);
    timer1->stop();

    m_nameLenProcessPtr->decompress(m_dnamelen.size, m_dnamelen.pdata);

    timer1->start("NameTime");
    m_nameProcessPtr->decompress(m_dname.size, m_dname.pdata);
    timer1->stop();

    if (m_paramPtr->m_alignType) 
    { // use reference
        if (m_paramPtr->m_dependent == DEPENDENT_SEQONQUAL) {
            timer1->start("QualTime");
            m_qualProcessPtr->decompress(m_dqual.size, m_dqual.pdata);
            timer1->stop();
            timer1->start("SeqTime");
            m_seqProcessPtr->decompress(m_dseq.size, m_dseq.pdata);
            timer1->stop();
        } else {
            timer1->start("SeqTime");
            m_seqProcessPtr->decompress(m_dseq.size, m_dseq.pdata);
            timer1->stop();
            timer1->start("QualTime");
            m_qualProcessPtr->decompress(m_dqual.size, m_dqual.pdata);
            timer1->stop();
        }
    } 
    else 
    {  // reference free
        if (m_paramPtr->m_seq_compress_mode == SEQ_REORDER && m_paramPtr->m_arcType == ARCTYPE_PE) 
        {
            timer1->start("SeqTime");
            m_refFreeSeqProcessPtr->decompress(m_dseq.size, m_dseq.pdata);
            timer1->stop();
            timer1->start("QualTime");
            m_qualProcessPtr->decompress(m_dqual.size, m_dqual.pdata);
            timer1->stop();
            pe_reorder_process();
        }
        else
        {
            timer1->start("CombineTime");
            m_combineProcessPtr->decompress(m_combine.size, m_combine.pdata);
            timer1->stop();
        }


        /*timer1->start("SeqTime");
        m_refFreeSeqProcessPtr->decompress(m_dseq.size, m_dseq.pdata);
        timer1->stop();
        timer1->start("QualTime");
        m_qualProcessPtr->decompress(m_dqual.size, m_dqual.pdata);
        timer1->stop();
        if (m_paramPtr->m_seq_compress_mode == SEQ_REORDER && m_paramPtr->m_arcType == ARCTYPE_PE) 
        {
            pe_reorder_process();
        }**/
    }
}

int DecodeWorker::copyReadDataNoComment(int i, char *strout) //复制数据用于输出
{
    char *pout = strout;
    *pout++ = '@';
    memcpy(pout, m_nameElPtr->m_bufPtr, m_namelenElPtr->m_LenArry[i]); //name
    pout += m_namelenElPtr->m_LenArry[i];
    m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
    *pout++ = '\n';

    memcpy(pout, m_seqElPtr->m_bufPtr, m_seqLenElPtr->m_LenArry[i]); //seq
    pout += m_seqLenElPtr->m_LenArry[i];
    m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];

    memcpy(pout, "\n+\n", 3); //comment
    pout += 3;

    memcpy(pout, m_qualElPtr->m_bufPtr, m_seqLenElPtr->m_LenArry[i]); //qual
    pout += m_seqLenElPtr->m_LenArry[i];
    m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
    *pout++ = '\n';

    return pout-strout;
}

int DecodeWorker::copyReadDataWithComment(int i, char *strout)
{
    char *pout = strout;
    *pout++ = '@';
    memcpy(pout, m_nameElPtr->m_bufPtr, m_namelenElPtr->m_LenArry[i]); //name
    pout += m_namelenElPtr->m_LenArry[i];
    *pout++ = '\n';

    memcpy(pout, m_seqElPtr->m_bufPtr, m_seqLenElPtr->m_LenArry[i]); //seq
    pout += m_seqLenElPtr->m_LenArry[i];
    m_seqElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
    *pout++ = '\n';

    *pout++ = '+';
    memcpy(pout, m_nameElPtr->m_bufPtr, m_namelenElPtr->m_LenArry[i]); //comment
    pout += m_namelenElPtr->m_LenArry[i];
    m_nameElPtr->m_bufPtr += m_namelenElPtr->m_LenArry[i];
    *pout++ = '\n';

    memcpy(pout, m_qualElPtr->m_bufPtr, m_seqLenElPtr->m_LenArry[i]); //qual
    pout += m_seqLenElPtr->m_LenArry[i];
    m_qualElPtr->m_bufPtr += m_seqLenElPtr->m_LenArry[i];
    *pout++ = '\n';

    return pout-strout;
}
