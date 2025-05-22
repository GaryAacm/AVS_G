#ifndef AVS_API_COMMANDPARSE_H
#define AVS_API_COMMANDPARSE_H

#include "Param.h"
#include <map>


class CommandParse
{
public:
    CommandParse();
    ~CommandParse();
    int parseOptFromCmd(int argc, char **argv);

private:
    void usage();
    void doCheckSetIndexOpt(int argc, char **argv);
    bool doCheckRefFile(const char *path);
    bool doCheckSrcFile(const char *path);
    void doMultipleCheck();
    void doCheckSetEncodeOpt(int argc, char **argv);
    void getOutFilePath();
    void doCheckSetDecodeOpt(int argc, char **argv);
    void setDecodeFilePath();
    bool isValidPath(string &str);
    void doSetExtractOpt(int argc, char **argv);
    void createOutFile();
private:
    bool m_force = false; //解压是否覆盖原文件
    bool m_bOutOriginPath = false; //是否保持输入文件同路径
    std::string m_strPrefix; //设定输出文件前缀
    Param *m_paramPtr = nullptr;
    std::string m_str_Multiple; //保存传入的多个文件路径，以','作为间隔
    std::string m_strExtract; //保存筛选参数
};




#endif //AVS_API_COMMANDPARSE_H