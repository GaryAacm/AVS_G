#include "AlignInfo.h"
#include "Param.h"

AlignInfo::AlignInfo()
{
    m_paramPtr = Param::GetInstance();
    m_cigarL.mem_malloc(128);
    m_cigarV.mem_malloc(128);
    m_seqarry.mem_malloc(m_maxlen);
    if(m_paramPtr->m_alignType == ALIGNTYPE_HASH)
    {
        int intlen = (m_maxlen>>4) + 1;
        m_seqbase.mem_malloc(intlen);
        m_revseqbase.mem_malloc(intlen);
        m_fqseed.mem_malloc(m_maxlen);
        m_revseqarry.mem_malloc(m_maxlen);
        m_revfqseed.mem_malloc(m_maxlen);
    }
}

AlignInfo::~AlignInfo()
{
}

void AlignInfo::init(int seqlen)
{
    if(m_maxlen < seqlen)
    {
        m_maxlen = seqlen;
        m_seqarry.mem_realloc(seqlen);
        if(m_paramPtr->m_alignType == ALIGNTYPE_HASH)
        {
            int intlen = (seqlen>>4) + 1;
            m_seqbase.mem_realloc(intlen);
            m_revseqbase.mem_realloc(intlen);
            m_fqseed.mem_realloc(seqlen);
            m_revseqarry.mem_realloc(seqlen);
            m_revfqseed.mem_realloc(seqlen);
        }
    }
    
    clear();
}

void AlignInfo::clear()
{
    m_bRev = false;
    m_misnum = 0xffffffff;
    m_seqlen = 0;
    h_degelen = 0;
    m_seedextendcount = 0;
    m_alignpos = 0;
    m_cigarL.clear();
    m_cigarV.clear();
    m_seqbase.clear();
    m_fqseed.clear();
    m_revseqarry.clear();
    m_revseqbase.clear();
    m_revfqseed.clear();
    m_seqarry.clear();
    m_refarry.clear();
    memset(m_hashseed, 0, sizeof(hashseed)* MAX_ALIGN_SEED);
}