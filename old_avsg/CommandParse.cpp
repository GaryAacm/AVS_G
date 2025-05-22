#include "CommandParse.h"
#include "ArcFile.h"
#include <fcntl.h>

CommandParse::CommandParse()
{
    m_paramPtr = Param::GetInstance();
}

CommandParse::~CommandParse()
{
}

int CommandParse::parseOptFromCmd(int argc, char **argv)
{
    int opt = 0;
    while ((opt = getopt(argc, argv, "hicdt:-fP:1:2:o:pn:m:B:x:D:r:k:")) != -1)
    {
        switch (opt) 
        {
            case 'h':
                usage();
                break;
            case 'i':
                m_paramPtr->m_actionType = ACTIONTYPE_DOINDEX;
                break;
            case 'c':
                m_paramPtr->m_actionType = ACTIONTYPE_DOENCODE;
                break;
            case 'd':
                m_paramPtr->m_actionType = ACTIONTYPE_DODECODE;
                break;
            case 't':
                m_paramPtr->m_iThreadnum = atoi(optarg);
                if (m_paramPtr->m_iThreadnum > MAX_THREAD_NUM) 
                {
                    arc_stdout(level::WARN, "Thread num can not be more than 100");
                    m_paramPtr->m_iThreadnum = MAX_THREAD_NUM;
                }
                if(m_paramPtr->m_iThreadnum == 0)
                {
                    m_paramPtr->m_iThreadnum = 1;
                }
                break;
            case '-':
                m_paramPtr->m_infiletype = INFILETYPE_PIPE;
                break;
            case 'f':
                m_force = true;
                break;
            case 'P':
                m_paramPtr->m_outfiletype = OUTFILETYPE_PIPE;
                m_paramPtr->m_PipeOut = (PipeOut)atoi(optarg);
                break;
            case '1':
                m_paramPtr->m_arcType = ARCTYPE_SE;
                m_paramPtr->m_vecFqFile.emplace_back(optarg);
                break;
            case '2':
                m_paramPtr->m_arcType = ARCTYPE_PE;
                m_paramPtr->m_vecFqFile.emplace_back(optarg);
                break;
            case 'o':
                m_strPrefix.append(optarg);
                break;
            case 'p':
                m_bOutOriginPath = true;
                break;
            case 'n':
                m_paramPtr->m_checktype = (DataCheck)atoi(optarg);
                break;
            case 'm':
                m_paramPtr->m_arcType = ARCTYPE_MULTI;
                m_str_Multiple.append(optarg);
                break;
            case 'B':
            {
                int size = atoi(optarg);
                if(size)
                {
                    m_paramPtr->m_iBlockSize = size;
                }
            }
                break;
            case 'x':
                m_strExtract.append(optarg);
                if (m_strExtract.empty())
                {
                    usage();
                }
                m_paramPtr->m_actionType = ACTIONTYPE_EXTRACT;
                break;
            case 'D':
                m_paramPtr->m_dependent = (Dependent)atoi(optarg);
                break;
            case 'r':
            {
                int is_reorder = atoi(optarg);
                if (is_reorder)
                {
                    m_paramPtr->m_seq_compress_mode = SEQ_REORDER;
                }
            }
                break;
            case 'k':
            {
                int kmer_size = atoi(optarg);
                if (kmer_size)  // TODO check kmer size range
                {
                    m_paramPtr->m_kmer_size = kmer_size;
                }
            }
                break;
            default:
                usage();
        }
    }

    m_paramPtr->m_iBlockSize *= 1024*1024;

    auto logger = Profile::getInstance().getModule("main")->getLogger();
    std::string name, stmp;


    switch (m_paramPtr->m_actionType)
    {
    case ACTIONTYPE_DOINDEX:
        doCheckSetIndexOpt(argc, argv);
        break;
    case ACTIONTYPE_DOENCODE:
        doCheckSetEncodeOpt(argc, argv);
        stmp.append(m_paramPtr->m_vecFqFile[0]);
        getFileName(stmp);
        name.append("Encode_").append(stmp);
        break;
    case ACTIONTYPE_DODECODE:
        doCheckSetDecodeOpt(argc, argv);
        stmp.append(m_paramPtr->m_strArcFile);
        getFileName(stmp);
        name.append("Decode_").append(stmp);
        break;
    case ACTIONTYPE_EXTRACT:
        doSetExtractOpt(argc, argv);
        stmp.append(m_paramPtr->m_strArcFile);
        getFileName(stmp);
        name.append("Extract_").append(stmp);
        break;
    default:
        arc_stdout(level::ERROR, "Please Check Input Param");
        usage();
        break;
    }

    logger->setLogName(name);
    return 0;
}

void CommandParse::usage()
{
    printf("AVS_api v%d.%d.%d. \n", MAJOR_VERS, MINOR_VERS, STAGE_VERS);
    printf("The entropy coder is derived from Fqzcomp. The aligner is derived from BWA.\n");
    printf("FASTA index is optional. Only Fqzcomp will be called if no index given or align ratio too low.\n");
    printf("For PE data, please compress each pair of files at once.\n\n");

    printf("To build HASH index:\n\tAVS_api -i <ref.fa>\n");
    printf("OR to build BWA index:\n\tAVS_api -q -i <ref.fa>\n");
    printf("To compress:\n\tAVS_api -c [options] [ref.fa] -1 <input_file> [-2 <input_file2>] -o <compress_prefix>\n");
    printf("To decompress:\n\tAVS_api -d <***.arc> -o <fastq_prefix> [options]\n\n");
    printf("\t          -t INT       Thread num for multi-threading, default as [1]\n");
    printf("\t          -h           Get instructions of this software\n");
    printf("\t          -f           Force overwrite target file\n");
    printf("\t          -P INT       Pipe out decompression result\n");
    printf("\t             1         Pipe out SE reads\n");
    printf("\t             2         Pipe out PE1 reads\n");
    printf("\t             3         Pipe out PE2 reads\n");
    printf("\t             4         Pipe out each pair of PE reads in order\n");
    printf("\t          -1           Set source file one\n");
    printf("\t          -2           Set source file two\n");
    printf("\t          -r INT       Set reference free compression mode, default as [0]\n");
    printf("\t            0          reference free compression preserve sequence order\n");
    printf("\t            1          reference free compression rearrange the sequence\n");
    printf("\t          -k INT       Set kmer size, default as [13]\n");
    printf("\t          -o           Set out file name\n");
    printf("\t          -p           Set out path to source file path\n");
    printf("\t          -n           Set data check type \n");
    printf("\t            0          set no check\n");
    printf("\t            1          set md5 check\n");
    printf("\t            2          set xxhash check\n");
    printf("\t          -m           Compress multi-file, use ',' to split file paths \n");
    printf("\t          -B           Set block size(MB), default as [50]\n");
    printf("\t          -x  string   extract reads\n");
    printf("\t              M:pattern   extract reads\n");
    printf("\t              H:1000   extract reads\n");
    printf("\t              S:1000   extract reads\n");
    printf("\t              L    extract reads\n");
    printf("\t              G:1,2   extract reads\n");
    printf("\t         -D INT        Set the dependency about seq and qual, default as [0]\n");
    printf("\t            0          seq and qual are independent\n");
    printf("\t            1          seq depend on qual\n");
    printf("\t            2          qual depend on seq\n");
    exit(0);
}

void CommandParse::doCheckSetIndexOpt(int argc, char **argv)
{
    int num = argc - optind;
    if(num != 1)
    {
        arc_stdout(level::ERROR, "Please Check Input Param");
        usage();
    }

    if(doCheckRefFile(argv[optind]))
    {
        m_paramPtr->m_strRefFile.append(argv[optind]);
    }
}

bool CommandParse::doCheckRefFile(const char *path)
{
    if(isValidRef(path))
    {
        if(isFileExist(path))
        {
            return true;
        }
        else
        {
            arc_stdout(level::ERROR, "The Ref file %s may be not exist", path);
        }
    }
    else
    {
        arc_stdout(level::ERROR, "The Ref file %s may be wrong type", path);
    }
    return 0;
}

bool CommandParse::doCheckSrcFile(const char *path)
{
    if(isValidSrc(path))
    {
        if(isFileExist(path))
        {
            return true;
        }
        else
        {
            arc_stdout(level::ERROR, "The Src file %s may be not exist", path);
        }
    }
    else
    {
        arc_stdout(level::ERROR, "The Src file %s may be wrong type", path);
    }
    return 0;
}

void CommandParse::doMultipleCheck()
{
    std::vector<std::string> vecfile;
    int pos = m_str_Multiple.find(',');
    if(pos != string::npos)
    {
        string strdeli(",");
        strSplit(m_str_Multiple, strdeli, vecfile);
    }
    else
    {
        fstream fin(m_str_Multiple.c_str(), std::ios::in);
        char tmpbuf[PATH_LEN]={0};
        while (!fin.eof())
        {
            fin.getline(tmpbuf, PATH_LEN);
            vecfile.emplace_back(tmpbuf);
        }
        if(vecfile[vecfile.size()-1].empty())
        {
            vecfile.pop_back(); //删除最后一行空行
        }
        fin.close();
    }

    for(string &str : vecfile)
    {
        doCheckSrcFile(str.c_str());
        m_paramPtr->m_vecFqFile.push_back(str);
    }
}

void CommandParse::getOutFilePath()
{   
    std::string strpath, strdir, strname;
    splitPath(m_paramPtr->m_vecFqFile[0].c_str(), strdir, strname);
    if(m_bOutOriginPath) //压缩路径和源文件同路径
    {
        strpath.append(strdir);
    }

    if(m_strPrefix.empty())
    {
        switch (m_paramPtr->m_arcType)
        {
        case ARCTYPE_SE:
            strname.append(".arc");
            break;
        case ARCTYPE_PE:
            strname.append(".PE.arc");
            break;
        case ARCTYPE_MULTI:
            strname.append(".Mul.arc");
            break;
        default:
            break;
        }
        m_strPrefix.append(strname);
    }
    else
    {
        splitPath(m_strPrefix.c_str(), strdir, strname);
        if(!strdir.empty())
        {
            if(m_bOutOriginPath) 
            {
                arc_stdout(level::ERROR, "the out path %s conflict with %s", m_strPrefix.c_str(), strpath.c_str());    
            }
            else
            {
                if(!isFileExist(strdir.c_str()))
                {
                    arc_stdout(level::ERROR, "the out path %s is not exist", strdir.c_str());
                }
            }
        }
        
        m_strPrefix.append(".arc");
    }

    strpath.append(m_strPrefix);
    m_paramPtr->m_strArcFile.append(strpath);
}

void CommandParse::doCheckSetEncodeOpt(int argc, char **argv)
{
    int num = argc - optind;
    if(num >1 || m_paramPtr->m_vecFqFile.empty())
    {
        arc_stdout(level::ERROR, "%s","Please Check Input Param");    
        usage();
    }

    if(num)
    {
        if(doCheckRefFile(argv[optind]))
        {
            m_paramPtr->m_strRefFile.append(argv[optind]);
        }
    }
    else
    {
        m_paramPtr->m_alignType = ALIGNTYPE_REFFREE;
    }

    if(m_paramPtr->m_arcType == ARCTYPE_MULTI)
    {
        doMultipleCheck();
    }
    else
    {
        for(string &str : m_paramPtr->m_vecFqFile)
        {
            doCheckSrcFile(str.c_str());
            if(str.find("gz") != string::npos)
            {
                m_paramPtr->m_totalfilesize += getFileSize(str.c_str())*5;
            }
            else
            {
                m_paramPtr->m_totalfilesize += getFileSize(str.c_str());
            }
        }
    }
    
    getOutFilePath();

    m_paramPtr->m_memBufPoolPtr = new MemBufPool(m_paramPtr->m_iThreadnum*2);
}

void CommandParse::doCheckSetDecodeOpt(int argc, char **argv)
{
    int num = argc - optind;
    if(num < 1 ||num >2)
    {
        arc_stdout(level::ERROR, "Please Check Input Param");
        usage();
    }

    if(num == 2)
    {
        if(doCheckRefFile(argv[optind]))
        {
            m_paramPtr->m_strRefFile.append(argv[optind]);
        }
        m_paramPtr->m_strArcFile.append(argv[optind+1]);
    }
    else if (num == 1)
    {
        m_paramPtr->m_strArcFile.append(argv[optind]);
    }

    if(!isFileExist(m_paramPtr->m_strArcFile.c_str()))
    {
        arc_stdout(level::ERROR, "The file %s may be not exist", m_paramPtr->m_strArcFile.c_str());
    }

    ArcFile arcfile;
    arcfile.readTail();
    if(m_paramPtr->m_PipeOut)
    {
        m_paramPtr->m_memBufPoolPtr = new MemBufPool(m_paramPtr->m_iThreadnum*2);
        if(m_paramPtr->m_arcType == ARCTYPE_SE)
        {
            m_paramPtr->m_PipeOut = PIPEOUT_SE;
        }
        else if(m_paramPtr->m_arcType == ARCTYPE_PE)
        {
            if(m_paramPtr->m_PipeOut <2 || m_paramPtr->m_PipeOut>4)
            {
                arc_stdout(level::ERROR, "the arc file is PE, please check Pipe param");
            }
        }
    }

    setDecodeFilePath();
    createOutFile();
}

bool CommandParse::isValidPath(string &str)
{
    if(!m_force)
    {
        if(isFileExist(str.c_str()))
        {
            arc_stdout(level::ERROR, "the out path %s is exist", str.c_str());
        }
    }
    return 0;
}

void CommandParse::setDecodeFilePath()
{
    std::string strpath, strdir, strname;
    splitPath(m_paramPtr->m_strArcFile.c_str(), strdir, strname);
    if(m_bOutOriginPath) //解压路径和压缩文件同路径
    {
        strpath.append(strdir);
    }
    if(!m_strPrefix.empty())
    {
        splitPath(m_strPrefix.c_str(), strdir, strname);
        if(!strdir.empty())
        {
            if(m_bOutOriginPath) //指定原路径
            {
                arc_stdout(level::ERROR, "the out path %s conflict with %s", strdir.c_str(), strpath.c_str());    
            }
            else
            {
                if(!isFileExist(strdir.c_str()))
                {
                    arc_stdout(level::ERROR, "the out path %s is not exist", strdir.c_str());
                }
            }
        }

        strpath.append(m_strPrefix);

        switch (m_paramPtr->m_arcType)
        {
        case ARCTYPE_PE:
        {
            m_paramPtr->m_vecFqFile[0] = strpath;
            m_paramPtr->m_vecFqFile[1] = strpath;
            m_paramPtr->m_vecFqFile[0].append("_1.fastq");
            m_paramPtr->m_vecFqFile[1].append("_2.fastq");
        }
            break;
        case ARCTYPE_SE:
            m_paramPtr->m_vecFqFile[0] = strpath;
            m_paramPtr->m_vecFqFile[0].append(".fastq");
            break;
        case ARCTYPE_MULTI:
        {
            auto itor = m_paramPtr->m_hashMultiple.begin();
            for(;itor!=m_paramPtr->m_hashMultiple.end();itor++)
            {
                int idx = *itor;
                m_paramPtr->m_vecFqFile[idx] = strpath;
                m_paramPtr->m_vecFqFile[idx].append(std::to_string(idx)).append(".fastq");
            }
            break;
        }
        default:
            break;
        }
    }
    else
    {
        if(!strpath.empty())
        {
            for(string &str : m_paramPtr->m_vecFqFile)
            {
                str = strpath.append(str);
            }
        }
    }
    

    for(string &str : m_paramPtr->m_vecFqFile)
    {
        isValidPath(str);
    }
}

void CommandParse::doSetExtractOpt(int argc, char **argv)
{
    int num = argc - optind;
    if(num < 1 ||num >2)
    {
        arc_stdout(level::ERROR, "Please Check Input Param");
        usage();
    }

    if(num == 2)
    {
        if(doCheckRefFile(argv[optind]))
        {
            m_paramPtr->m_strRefFile.append(argv[optind]);
        }
        m_paramPtr->m_strArcFile.append(argv[optind+1]);
    }
    else if (num == 1)
    {
        m_paramPtr->m_strArcFile.append(argv[optind]);
    }

    if(!isFileExist(m_paramPtr->m_strArcFile.c_str()))
    {
        arc_stdout(level::ERROR, "The file %s may be not exist", m_paramPtr->m_strArcFile.c_str());
    }

    ArcFile arcfile;

    switch (m_strExtract[0])
    {
    case 'M': //正则筛选
        arcfile.readTail();
        m_paramPtr->m_strPattern = m_strExtract.substr(2);
        m_paramPtr->m_extracttype = EXTRACT_MATCH;
        break;
    case 'H': //提取前n条read
    {
        std::map<uint32_t, BlockParamInfo> mapBParam;
        arcfile.readBlockParam(mapBParam); //遍历文件获取block信息
        arcfile.readTail();
        int totalnum = std::stoi(m_strExtract.substr(2));
        m_paramPtr->m_extracttype = EXTRACT_HEAD;
        int i = 0, readcnt = 0;
        std::vector<BlockMetaInfo> vectemp;
        while (totalnum>0)
        {
            readcnt = mapBParam[i].readcount;
            if(totalnum > readcnt)
            {
                m_paramPtr->m_vecReadCnt.push_back(readcnt);
                totalnum -= readcnt;
            }
            else
            {
                m_paramPtr->m_vecReadCnt.push_back(totalnum);
                totalnum -= totalnum;
            }
            vectemp.emplace_back(m_paramPtr->m_vecBlockMeta[mapBParam[i].writeidx]);
            i++;
        }
        m_paramPtr->m_vecBlockMeta.swap(vectemp);

        if(m_paramPtr->m_iThreadnum > m_paramPtr->m_vecReadCnt.size())
        {
            m_paramPtr->m_iThreadnum = m_paramPtr->m_vecReadCnt.size();
        }
        break;
    }
    case 'S': //随机提取n条read
    {
        arcfile.readTail();
        int totalnum = std::stoi(m_strExtract.substr(2));
        int min = totalnum < m_paramPtr->m_iBlocknum ? totalnum : m_paramPtr->m_iBlocknum;
        int k = totalnum / min;
        int mod = totalnum % min;
        int i = 0;
        for(;i<min;i++)
        {
            if(mod)
            {
                m_paramPtr->m_vecReadCnt.push_back(k+1);
                mod--;
            }
            else
            {
                m_paramPtr->m_vecReadCnt.push_back(k);
            }
        }
        if(m_paramPtr->m_iThreadnum > m_paramPtr->m_vecReadCnt.size())
        {
            m_paramPtr->m_iThreadnum = m_paramPtr->m_vecReadCnt.size();
        }
        m_paramPtr->m_extracttype = EXTRACT_SRAND;
        break;
    }
    case 'L': //显示压缩的多文件序号
        m_paramPtr->m_extracttype = EXTRACT_LIST;
        break;
    case 'G': //提取相应序号的文件
    {
        m_paramPtr->m_extracttype = EXTRACT_GETFILES;
        vector<string> vec_idx;
        string strdeli(",");
        strSplit(m_strExtract.substr(2), strdeli, vec_idx);
        for(int i=0;i<vec_idx.size();i++)
        {
            int k = std::stoi(vec_idx[i]);
            m_paramPtr->m_hashMultiple.insert(k);
        } 
        arcfile.readTail();       
    }
        break;
    default:
        break;
    }

    m_paramPtr->m_memBufPoolPtr = new MemBufPool(m_paramPtr->m_iThreadnum*2);
    
    setDecodeFilePath();
    createOutFile();
}

void CommandParse::createOutFile()
{
    if(m_paramPtr->m_arcType == ARCTYPE_MULTI|| m_paramPtr->m_PipeOut)
    {
        return;
    }
    m_paramPtr->m_outfd[0] = open(m_paramPtr->m_vecFqFile[0].c_str(), O_RDWR|O_CREAT, 0666);
    if(m_paramPtr->m_outfd[0] < 0)
    {
        arc_stdout(level::ERROR, "create %s fail", m_paramPtr->m_vecFqFile[0].c_str());
    }
    if(m_paramPtr->m_arcType == ARCTYPE_PE)
    {
        m_paramPtr->m_outfd[1] = open(m_paramPtr->m_vecFqFile[1].c_str(), O_RDWR|O_CREAT, 0666);
        if(m_paramPtr->m_outfd[1] < 0)
        {
            arc_stdout(level::ERROR, "create %s fail", m_paramPtr->m_vecFqFile[1].c_str());
        }
    }
}
