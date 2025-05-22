#ifndef AVS_API_REFFREESEQPROCESS_H
#define AVS_API_REFFREESEQPROCESS_H

#include "IProcess.h"
#include "Encap.h"
#include "IRef.h"
#include "graph.h"
// #include <string_view>
#include "string-view-lite/include/nonstd/string_view.hpp"
#include <unordered_map>
#include <numeric>
#include <queue>
#include <array>
#include <boost/filesystem.hpp>
#include <thread>
#include <random>

// 注意使用Element对象需要前向声明
class BlockParamElement;

class IWorker;

class SeqElement;

class SeqLenElement;

class QualEmement;

static constexpr size_t MAX_SHORT_READ_LEN = 511; // 短读序列长度最高为511

class RefFreeSeqProcess : public IProcess {
public:
    explicit RefFreeSeqProcess(IWorker *ptr);

    ~RefFreeSeqProcess() override;

    int compress(char *outptr) override;

    int decompress(int inlen, char *inptr) override;

private:
    int calcMd5(char *outptr);
    bool checkData();
    void clear_all(); // 在执行压缩/解压之前清除类中所有成员变量
    int do_compress_seq(char *outptr);
    int internal_compress(char *output);

    // another implement
    // void construct_reference(); // build index and construct zero_error ref
    // void alignment();           // index other reads and reference, align reads to ref (use copMEM)

    // kuafu-cpp implement
    void build_minimizer_index_and_alignment(Graph &graph);
    std::string get_seq(uint32_t read_id);
    void graph_encode(Graph& graph);
    int encode_output(char * output);

    // long read compress
    int long_read_compress(char * outptr);
    int long_read_compress_implement(char * outptr);

    void clean_working_directory();

    SeqElement *m_seqElPtr = nullptr;
    SeqLenElement *m_seqLenElPtr = nullptr;
    QualEmement *m_qualElPtr = nullptr;
    Param *m_param = nullptr;
    std::shared_ptr<Module> m_profile;
    boost::filesystem::path working_path;

    uint32_t max_seq_len = 0;
    bool is_long_read_compress = false;
    unsigned char m_szmd5[16];
    XXH128_hash_t m_hashval;
    std::vector<nonstd::string_view> seq_view;
    std::vector<std::string> N_reconstruct;
    std::vector<std::string> N_encode;
    std::vector<bool> is_contain_N;

public:
    std::vector<uint32_t> order; // order[旧序列顺序ID]储存新的序列顺序ID
    std::vector<uint8_t> pe_id;  // 重排序列的双端文件解压时需要pe_id区分不同文件的序列
};


#endif //AVS_API_REFFREESEQPROCESS_H