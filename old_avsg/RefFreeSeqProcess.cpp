#include "RefFreeSeqProcess.h"
#include "IWorker.h"
#include "dna.h"
#include "environment.h"
#include <cstdlib>
#include <omp.h>
#include "clr.h"
#include "simple_model.h"

#define DEFINE_PATH_MACRO                                                                                                        \
boost::filesystem::path                                                                                                          \
        extra_read_path = working_path / "extra_read.txt", extra_read_comp_path = working_path / "extra_read.bsc",               \
        read_path = working_path / "read.txt", read_comp_path = working_path / "read.bsc",                                       \
        overlap_dir_path = working_path / "overlap_dir.txt", overlap_dir_comp_path = working_path / "overlap_dir.bsc",           \
        strand_dir_path = working_path / "strand_dir.txt", strand_dir_comp_path = working_path / "strand_dir.bsc",               \
        mismatch_path = working_path / "mismatch.txt", mismatch_comp_path = working_path / "mismatch.bsc",                       \
        tree_path = working_path / "tree.txt", tree_comp_path = working_path / "tree.bsc",                                       \
        overlap_offset_path = working_path / "overlap_off.bin", overlap_offset_comp_path = working_path / "overlap_off.bin.bsc", \
        pe_id_path = working_path / "pe.txt", pe_id_comp_path = working_path / "pe.bsc",                                         \
        tree_read_id_path = working_path / "tree_read_id.bin", tree_read_id_comp_path = working_path / "tree_read_id.bin.bsc"    \


RefFreeSeqProcess::RefFreeSeqProcess(IWorker *ptr) 
{
    m_seqElPtr = ptr->getseqElPtr();
    m_seqLenElPtr = ptr->getseqLenElPtr();
    m_qualElPtr = ptr->getqualElPtr();
    m_param = Param::GetInstance();
    m_profile = ptr->m_profile;

    std::random_device device;
    std::mt19937_64 gen(device());
    std::uniform_int_distribution<uint32_t> dis(0, std::numeric_limits<uint32_t>::max());
    while (true) {
        working_path = boost::filesystem::path("tmp") / std::to_string(dis(gen));
        if (!boost::filesystem::exists(working_path)) {
            boost::filesystem::create_directories(working_path);
            break;
        }
    }
}

RefFreeSeqProcess::~RefFreeSeqProcess() {

}

std::string bsc_compress(const char *input, const char *output) {
    if (!boost::filesystem::exists(bsc_path)){
        throw std::runtime_error("BSC_PATH isn't a valid path");
    }
    boost::filesystem::path input_path(input);
    boost::filesystem::path output_path(output);
    std::string bsc_command = std::string(bsc_path) + " e "
                              + boost::filesystem::absolute(input_path).c_str() + " "
                              + boost::filesystem::absolute(output_path).c_str() +
                              " -p -e2 ";
    return bsc_command;
}

void bsc_decompress(const char *input, const char *output) {
    if (!boost::filesystem::exists(bsc_path)){
        throw std::runtime_error("BSC_PATH isn't a valid path");
    }
    boost::filesystem::path input_path(input);
    boost::filesystem::path output_path(output);
    std::string bsc_command = std::string(bsc_path) + " d "
                              + boost::filesystem::absolute(input_path).c_str() + " "
                              + boost::filesystem::absolute(output_path).c_str() + " "
                              + " > /dev/null"; // supress output
    int status = std::system(bsc_command.c_str());
    if (status != 0) throw std::runtime_error("Error occurred during bsc decompress.");
}

template<typename T>
static inline void write_data_to_binary_file(const T *data, size_t data_size, const char *filename) {
    std::ofstream file;
    file.open(filename, std::ios::out | std::ios::binary);
    if (data_size > 0) {
        file.write(reinterpret_cast<const char *>(data), sizeof(T) * data_size);
    }
}

template<typename T>
static inline void write_vector_to_binary_file(const std::vector<T> &data, const std::string &filename) {
    write_data_to_binary_file(data.data(), data.size(), filename.c_str());
}

template<typename T>
static inline void read_vector_from_binary_file(std::vector<T> &data, const std::string &filename) {
    std::ifstream file;
    file.open(filename, std::ios::in | std::ios::binary);
    file.seekg(0, std::ios::end);
    auto size = file.tellg();
    file.seekg(0);
    data.resize(size / sizeof(T));
    file.read(reinterpret_cast<char *>(data.data()), size);
}

/// 整数压缩器
/// \param data
/// \param compressed_output 提供的compress_output可以为空数组
/// \return
void Integer_compress(const std::vector<u32> &data, std::vector<char> &compressed_output) {
    compressed_output.resize(data.size() * sizeof(u32) + 1024);
    char * output = compressed_output.data();
    const u32 origin_data_size = data.size();
    std::memcpy(output, &origin_data_size, 4);
    output += 4;

    RangeCoder rc;
    rc.output(output);
    rc.StartEncode();

    SIMPLE_MODEL<256> len1;
    SIMPLE_MODEL<256> len2;
    SIMPLE_MODEL<256> len3;
    SIMPLE_MODEL<256> len4;
    SIMPLE_MODEL<2> same_len;
    u32 m_last = 0;

    for (uint32_t i = 0; i < data.size(); i++) {
        if (data[i] != m_last) {
            same_len.encodeSymbol(&rc, 0);
            len1.encodeSymbol(&rc, data[i] & 0xff);
            len2.encodeSymbol(&rc, (data[i] >> 8) & 0xff);
            len3.encodeSymbol(&rc, (data[i] >> 16) & 0xff);
            len4.encodeSymbol(&rc, (data[i] >> 24) & 0xff);
            m_last = data[i];
        } else {
            same_len.encodeSymbol(&rc, 1);
        }
    }

    rc.FinishEncode();
    compressed_output.resize(rc.size_out() + 4);
}

/// 整数解压器
/// \param compress_data
/// \param decompress_output 提供的decompress_output必须为空数组
void Integer_decompress(const std::vector<char> &compress_data, std::vector<u32> &decompress_output) {
      u32 origin_data_size = 0;
      char* input = const_cast<char*>(compress_data.data());
      std::memcpy(&origin_data_size, input, 4);
      input += 4;

      RangeCoder rc;
      rc.input(input);
      rc.StartDecode();

      SIMPLE_MODEL<256> len1;
      SIMPLE_MODEL<256> len2;
      SIMPLE_MODEL<256> len3;
      SIMPLE_MODEL<256> len4;
      SIMPLE_MODEL<2> same_len;
      u32 m_last = 0;

      decompress_output.reserve(origin_data_size);
      int l1=0,l2=0,l3=0,l4=0;
      for (int i = 0; i < origin_data_size; i++) {
          if (same_len.decodeSymbol(&rc)) {
              decompress_output.push_back(m_last);
          } else {
              l1 = len1.decodeSymbol(&rc);
              l2 = len2.decodeSymbol(&rc);
              l3 = len3.decodeSymbol(&rc);
              l4 = len4.decodeSymbol(&rc);
              m_last = l1 + (l2 << 8) + (l3<<16) + (l4<<24);
              decompress_output.push_back(m_last);
          }
      }

      rc.FinishDecode();
}

static constexpr char table[4] = {'A', 'C', 'G', 'T'};
struct Count {
    std::array<uint32_t, 4> count{};

    void add(char c) {
        switch (c) {
            case 'A':
                count[0]++;
                break;
            case 'C':
                count[1]++;
                break;
            case 'G':
                count[2]++;
                break;
            case 'T':
                count[3]++;
                break;
            default :
                break;
        }
    }

    [[nodiscard]] char get() const noexcept {
        return table[std::distance(count.begin(), std::max_element(count.begin(), count.end()))];
    }
};

struct Minimizer {
    u64 key;
    u64 read_id: 32, position: 31, strand_id: 1;

    Minimizer(u64 key, u32 read_id, u32 pos, u8 strand_id)
            : key(key), read_id(read_id), position(pos), strand_id(strand_id) {}

    Minimizer() : key(0), read_id(0), position(0), strand_id(0) {}

    friend bool operator<(const Minimizer &m1, const Minimizer &m2) {
        return std::tie(m1.key, m1.position, m1.read_id, m1.strand_id) <
               std::tie(m2.key, m2.position, m2.read_id, m2.strand_id);
    }

    friend bool operator==(const Minimizer &m1, const Minimizer &m2) {
        return m1.key == m2.key &&
               m1.read_id == m2.read_id &&
               m1.position == m2.position &&
               m1.strand_id == m2.strand_id;
    }
};

struct Match {
    u8 start_strand_id: 4, end_strand_id: 4;
    u8 mismatch;
    int16_t offset; // offset > 0 : start_suffix <-> end_prefix ; offset < 0 : start_prefix <-> end_suffix

    explicit Match() : start_strand_id(0), end_strand_id(0), mismatch(0), offset(0) {}

    Match(u32 start_strand_id, u32 end_strand_id, u32 mismatch, int offset) :
            start_strand_id(start_strand_id), end_strand_id(end_strand_id), mismatch(mismatch), offset(offset) {}
};

struct AlignmentElement {
    u32 target_read_id;
    Match match;

    AlignmentElement(u32 target_read_id, const Match &match) : target_read_id(target_read_id), match(match) {}

    friend bool operator<(const AlignmentElement &a, const AlignmentElement &b) {
        int a_off = std::abs(a.match.offset), b_off = std::abs(b.match.offset);
        u32 a_mismatch = a.match.mismatch, b_mismatch = b.match.mismatch;
        return std::tie(a_off, a_mismatch) < std::tie(b_off, b_mismatch);
    }

    friend bool operator>(const AlignmentElement &a, const AlignmentElement &b) {
        int a_off = std::abs(a.match.offset), b_off = std::abs(b.match.offset);
        u32 a_mismatch = a.match.mismatch, b_mismatch = b.match.mismatch;
        return std::tie(a_off, a_mismatch) > std::tie(b_off, b_mismatch);
    }
};

struct MinimizerHashTable {
    /// 仅记录每个position区间所在的位置(相邻区间之间重复的位置只记录一次),并且不需要记录position,必要的时候直接从sort_minimizer取出来即可.
    /// 即:
    ///   0 1 2 3   4 5 6   7 8   9 10   11
    /// [ o o o o | o o o | o o | o o ] end
    /// PositionMap = { 0, 4, 7, 9, 11 } => 四个区间: [0,4),[4,7),[7,9),[9,11)
    using PositionMap = std::vector<size_t>;

private:
    std::vector<Minimizer> sort_minimizer_array;
    std::unordered_map<u64, PositionMap> range_table; /// minimizer_key -> vector<pos -> range>
    std::vector<std::vector<u32>> minimizer_access_index; /// 用于快速访问一个序列所在的minimizer
public:
    explicit MinimizerHashTable(size_t read_count, size_t minimizer_num) {
        sort_minimizer_array.reserve(read_count * minimizer_num);
        minimizer_access_index.resize(read_count);
        for (auto& item: minimizer_access_index) item.reserve(minimizer_num);
    }

    inline void sort_minimizer() {
        sort_minimizer_array.shrink_to_fit();
        std::sort(sort_minimizer_array.begin(), sort_minimizer_array.end(), std::less<Minimizer>{});

        for (u32 access_index = 0; access_index < sort_minimizer_array.size(); ++access_index){
            minimizer_access_index[sort_minimizer_array[access_index].read_id].push_back(access_index);
        }

        u64 key = sort_minimizer_array[0].key;
        std::vector<size_t> key_range;
        for (size_t i = 1; i < sort_minimizer_array.size(); i++) {
            if (sort_minimizer_array[i].key != key) {
                key_range.push_back(i);
                key = sort_minimizer_array[i].key;
            }
        }
        key_range.push_back(sort_minimizer_array.size());

        PositionMap position_map;
        size_t begin = 0;
        for (size_t end:key_range) {
            position_map.clear();
            position_map.push_back(begin);
            u32 position = sort_minimizer_array[begin].position;
            for (size_t i = begin + 1; i < end; i++) {
                u32 temp = sort_minimizer_array[i].position;
                if (temp != position) {
                    position_map.push_back(i);
                    position = temp;
                }
            }
            position_map.push_back(end);
            u64 hash_key = sort_minimizer_array[begin].key;
            range_table.emplace(hash_key, position_map);
            begin = end;
        }
    }

    [[nodiscard]] inline const std::vector<u32>& get_minimizer_access_index(u32 read_id) const {
        return minimizer_access_index[read_id];
    }

    // 该函数不会检查是否存在key,使用前就应该确保key存在
    [[nodiscard]] inline const PositionMap &get(const Minimizer &m) const {
        return (*(range_table.find(m.key))).second;
    }

    [[nodiscard]] inline std::vector<Minimizer> &data() {
        return sort_minimizer_array;
    }

    [[nodiscard]] inline const std::vector<Minimizer> &data() const {
        return sort_minimizer_array;
    }
};

template<typename T1, typename T2>
inline void concat_vector(std::vector<T1> &target, const std::vector<T2> &source) {
    target.insert(target.end(), source.begin(), source.end());
}

template<typename IntegerType>
static inline u64 invertible_hash(IntegerType key) {
    IntegerType mask = std::numeric_limits<IntegerType>::max();
    key = (~key + (key << 21u)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24u;
    key = ((key + (key << 3u)) + (key << 8u)) & mask; // key * 265
    key = key ^ key >> 14u;
    key = ((key + (key << 2u)) + (key << 4u)) & mask; // key * 21
    key = key ^ key >> 28u;
    key = (key + (key << 31u)) & mask;
    return key;
}

std::vector<Minimizer> get_minimizer(size_t read_id, const nonstd::string_view &seq, size_t k, size_t w, size_t m) {
    constexpr u64 MAX_U64 = std::numeric_limits<u64>::max();
    std::vector<Minimizer> minimizer_set;
    std::vector<u64> kmers, rev_kmers;

    for (u64 i = 0; i < seq.size() - k + 1; i++) {
        u64 kmer = kmer_to_u64(seq.substr(i, k));
        u64 rev_kmer = get_reverse_kmer(kmer, k);
        kmers.push_back(invertible_hash(kmer));
        rev_kmers.push_back(invertible_hash(rev_kmer));
    }
    std::reverse(rev_kmers.begin(), rev_kmers.end());

    // 起始的w个连续的kmer和reverse_kmer
    u64 begin_pos = 0, end_pos = w;
    u64 minimizer_pos = 0;
    u64 minimizer_key = MAX_U64;
    auto compute_minimizer = [&]() -> void {
        minimizer_key = MAX_U64;
        for (u64 pos = begin_pos; pos < end_pos; ++pos) {
            auto temp = std::min(kmers[pos], rev_kmers[pos]);
            if (minimizer_key > temp) {
                minimizer_key = temp;
                minimizer_pos = pos;
            }
        }
        for (u64 pos = begin_pos; pos < end_pos; ++pos) {
            // 由于我们强制kmer_size为奇数,所以这里不存在kmer==rev_kmer的情况.
            if (kmers[pos] < rev_kmers[pos] && kmers[pos] == minimizer_key) {
                minimizer_set.emplace_back(minimizer_key, read_id, pos, 0);
                if (minimizer_set.size() >= m) return;
            } else if (kmers[pos] > rev_kmers[pos] && rev_kmers[pos] == minimizer_key) {
                minimizer_set.emplace_back(minimizer_key, read_id, pos, 1);
                if (minimizer_set.size() >= m) return;
            }
        }
    };
    while (end_pos <= kmers.size()) {
        if (begin_pos == 0) { // first window
            compute_minimizer();
        } else { // 通过上一次的结果加速
            if (minimizer_pos < begin_pos) { // 上一次的minimizer位置不在当前区域内
                compute_minimizer();
            } else {
                auto temp = std::min(kmers[end_pos - 1], rev_kmers[end_pos - 1]);
                if (temp <= minimizer_key) { // 如果temp>minimizer_key,那么当前区间内的minimizer依然不变,不需要更新
                    minimizer_key = temp;
                    minimizer_pos = end_pos - 1;
                    if (kmers[end_pos - 1] < rev_kmers[end_pos - 1]) {
                        minimizer_set.emplace_back(minimizer_key, read_id, end_pos - 1, 0);
                    } else {
                        minimizer_set.emplace_back(minimizer_key, read_id, end_pos - 1, 1);
                    }
                }
            }
        }
        if (minimizer_set.size() >= m) return minimizer_set;
        begin_pos++;
        end_pos++;
    }
    return minimizer_set;
}

/// 对于 window_size <= 0 的, 直接提取 kmer 即可
std::vector<Minimizer> get_minimizer2(size_t read_id, const nonstd::string_view &seq, size_t k, size_t m){
    std::vector<Minimizer> minimizer_set;
    size_t kmer_num = seq.size() - k + 1; 
    size_t begin_pos = 0, end_pos = 0;
    if (kmer_num <= m) {
        begin_pos = 0; end_pos = kmer_num;
    } else {
        begin_pos = kmer_num / 2 - m / 2; end_pos = begin_pos + m ;
    }
    for (size_t pos = begin_pos; pos < end_pos; ++pos) {
        u64 kmer = kmer_to_u64(seq.substr(pos, k));
        u64 rev_kmer = get_reverse_kmer(kmer, k);
        kmer = invertible_hash(kmer);
        rev_kmer = invertible_hash(rev_kmer);
        if (kmer < rev_kmer) {
            minimizer_set.emplace_back(kmer, read_id, pos, 0);
        } else {
            minimizer_set.emplace_back(rev_kmer, read_id, seq.size() - k - pos , 1);
        }
    }
    return minimizer_set;
}

static inline void reverse_complement(std::string &str) {
    std::reverse(str.begin(), str.end());
    for (auto &ch:str) {
        switch (ch) {
            case 'A':
                ch = 'T';
                break;
            case 'T':
                ch = 'A';
                break;
            case 'C':
                ch = 'G';
                break;
            case 'G':
                ch = 'C';
                break;
            default:
                break;
        }
    }
}

static inline void reverse_complement(char *data, int len) {
    for (int i = 0; i < len; ++i) {
        switch (data[i]) {
            case 'A':
                data[i] = 'T';
                break;
            case 'T':
                data[i] = 'A';
                break;
            case 'C':
                data[i] = 'G';
                break;
            case 'G':
                data[i] = 'C';
                break;
            default:
                break;
        }
    }
    for (int i = 0; i < len / 2; ++i) {
        std::swap(data[i], data[len - 1 - i]);
    }
}

void RefFreeSeqProcess::build_minimizer_index_and_alignment(Graph &graph) {
    size_t minimizer_num = 4; // 10;
    size_t kmer_size = m_param->m_kmer_size;
    size_t alignment_result_max = 2; // 5;
    size_t max_match_times = 5;      // 1000;
    double min_se_overlap_rate = 0.5;  // 0.3;

    auto alignment_two_read = [&](const Minimizer &u, const Minimizer &v, size_t max_se_mismatch) -> Match {
        Match match{};
        match.offset = (int) (u.position) - (int) (v.position);
        match.start_strand_id = u.strand_id;
        match.end_strand_id = v.strand_id;
        u32 mismatch = 0;
        std::string u_seq, v_seq;
        u_seq.reserve(seq_view[u.read_id].size());
        v_seq.reserve(seq_view[v.read_id].size());
        for (auto ch: seq_view[u.read_id]) u_seq += ch;
        for (auto ch: seq_view[v.read_id]) v_seq += ch;
        if (u.strand_id == 1) reverse_complement(u_seq);
        if (v.strand_id == 1) reverse_complement(v_seq);
        if (match.offset >= 0) { // u_suffix - v.prefix
            size_t overlap = std::min(u_seq.size() - match.offset, v_seq.size());
            for (size_t pos = 0; pos < overlap; pos++) {
                if (u_seq[pos + match.offset] != v_seq[pos]) {
                    mismatch++;
                    if (mismatch > (u32) max_se_mismatch) break;
                }
            }
        } else {           // u_prefix - v.suffix
            size_t overlap = std::min(v_seq.size() + match.offset, u_seq.size());
            for (size_t pos = 0; pos < overlap; pos++) {
                if (u_seq[pos] != v_seq[pos - match.offset]) {
                    mismatch++;
                    if (mismatch > (u32) max_se_mismatch) break;
                }
            }
        }
        match.mismatch = mismatch;
        return match;
    };

    // index
    MinimizerHashTable index_result(seq_view.size(), minimizer_num);
#pragma omp parallel for schedule(dynamic, 10000)
    for (uint32_t seq_id = 0; seq_id < seq_view.size(); ++seq_id) {
        if (seq_view[seq_id].size() < kmer_size) continue;
        std::vector<Minimizer> ret;
        int window_size = static_cast<int>(seq_view[seq_id].size() / 2) - static_cast<int>(kmer_size);
        if (window_size <= 0) {
            ret = get_minimizer2(seq_id, seq_view[seq_id], kmer_size, minimizer_num);
        } else {
            // size_t window_size = seq_view[seq_id].size() / 2 - kmer_size;
            ret = get_minimizer(seq_id, seq_view[seq_id], kmer_size, window_size, minimizer_num);
        }
#pragma omp critical
        concat_vector(index_result.data(), ret);
    }
    index_result.sort_minimizer();
    const auto &data = index_result.data();

    // alignment
#pragma omp parallel for schedule(dynamic, 10000)
    for (u32 read_id = 0; read_id < seq_view.size(); ++read_id) {
        double error_rate = 0.04;
        size_t read_length = seq_view[read_id].size();
        size_t max_se_mismatch = read_length * error_rate;
        if (read_length < 100) max_se_mismatch = 4;
        if (read_length < kmer_size) continue;
        size_t min_se_overlap = read_length * min_se_overlap_rate;
        std::vector<AlignmentElement> buf;
        buf.reserve(alignment_result_max);
        const auto& M = index_result.get_minimizer_access_index(read_id);
        int alignment_times = 2;
        while (alignment_times --) {
            for (const auto &m_id : M) {
                const auto &m = data[m_id];
                const auto &range = index_result.get(m);
                std::vector<size_t> permutation(range.size() - 1);
                std::iota(permutation.begin(), permutation.end(), 0);
                std::sort(permutation.begin(), permutation.end(), [&](size_t a, size_t b) {
                    int a_pos = data[range[a]].position;
                    int b_pos = data[range[b]].position;
                    return std::abs((int) (m.position) - a_pos) <
                           std::abs((int) (m.position) - b_pos);
                });

                u32 emplace_op_time = 0;
                for (u32 index:permutation) {
                    auto begin = range[index], end = range[index + 1];
                    int offset = (int) (m.position) - (int) (data[begin].position);
                    if ((read_length - std::abs(offset)) < min_se_overlap) break;

                    auto data_begin = data.begin() + begin, data_end = data.begin() + end;
                    auto itr = std::upper_bound(data_begin, data_end, read_id,
                                                [](u32 read_id, const Minimizer &m) {
                                                    return read_id < m.read_id;
                                                });

                    u32 match_calculation_time = 0;

                    while (itr != data_end) {
                        const auto &n = *itr;
                        auto ret = std::find_if(buf.begin(), buf.end(), [&](const AlignmentElement &e) {
                            return e.target_read_id == n.read_id;
                        });
                        if (ret != buf.end()) {
                            if (std::abs(offset) < std::abs(ret->match.offset)) {
                                auto match = alignment_two_read(m, n, max_se_mismatch);
                                if ( match.mismatch <= max_se_mismatch) {
                                    ret->match = match;
                                }
                                match_calculation_time++;
                            }
                        } else {
                            auto match = alignment_two_read(m, n, max_se_mismatch);
                            if (match.mismatch <= max_se_mismatch) {
                                AlignmentElement ele(n.read_id, match);
                                if (buf.size() < alignment_result_max) {
                                    buf.push_back(ele);
                                } else {
                                    auto max_itr = std::max_element(buf.begin(), buf.end());
                                    if (ele < *max_itr) *max_itr = ele;
                                }
                                emplace_op_time++;
                                if (emplace_op_time >= alignment_result_max) goto alignment_end;
                            }
                            match_calculation_time++;
                        }
                        if (match_calculation_time >= max_match_times) break;
                        itr++;
                    }
                }
                alignment_end:;
            }
            if (buf.size() >= alignment_result_max) break;
            error_rate += 0.04;
            if (read_length >= 100) max_se_mismatch = read_length * error_rate;
            else {
                max_se_mismatch += 4;
                if (max_se_mismatch > (read_length / 2)) break;
            }
        }

        auto size = std::min(buf.size(), alignment_result_max);
#pragma omp critical
        for (u32 i = 0; i < size; i++) {
            u32 target_read_id = buf[i].target_read_id;
            const auto &m = buf[i].match;
            Graph::G::Arc e;
            if (m.offset >= 0) {
                e = graph.g.addArc(Graph::G::nodeFromId(read_id), Graph::G::nodeFromId(target_read_id));
                graph.strand_map[e] = Strand(m.start_strand_id, m.end_strand_id);
            } else {
                e = graph.g.addArc(Graph::G::nodeFromId(target_read_id), Graph::G::nodeFromId(read_id));
                graph.strand_map[e] = Strand(m.end_strand_id, m.start_strand_id);
            }
            graph.weight_map[e] = Weight(std::abs(m.offset), m.mismatch);
        }
    }

}

/// child_seq 和 parent_seq 已经按照其方向纠正.
/// 对于后缀-前缀重叠,可能出现child_seq包含在parent_seq的情况,如下:
///     --------------------- parent_seq
///          ------------ child_seq. (虽然也是后缀-前缀重叠,但是child_seq和原序列长度不一致就会出现这种情况).
/// 此时不需要额外记录重叠以外的部分,当然偏移量也是要记录的.
/// 对于前缀-后缀重叠,可能出现parent_seq包含在child_seq的情况,如下:
///        ----------------       parent_seq
///   --------------------------   child_seq
////  <-A->                <-B->
/// 此时需要把A与B两部分按顺序保存下来,同时记录偏移量,解压的时候先根据偏移量解析A,然后利用
/// parent_seq和child_seq顺序比对解压完中间部分,最后利用已知的child_seq长度减去已解压的长度,再往下解压B部分即可.
std::string compute_overlap(const nonstd::string_view &child_seq, const nonstd::string_view &parent_seq,
                            u8 overlap_dir, int off) {
    std::string dna_string;
    if (overlap_dir == 0) {
        int leave_size = off + (int) child_seq.size() - (int) parent_seq.size();
        if (leave_size > 0) {
            dna_string = child_seq.substr(child_seq.size() - leave_size, leave_size).to_string();
        }
    } else {
        dna_string = child_seq.substr(0, off).to_string();
        int leave_size = (int) child_seq.size() - (off + (int) parent_seq.size());
        if (leave_size > 0) {
            dna_string += child_seq.substr(child_seq.size() - leave_size, leave_size).to_string();
        }
    }
    return dna_string;
}

inline std::string compute_mismatch_string(
        const nonstd::string_view &correct_seq, const nonstd::string_view &mismatch_seq, int offset) {
    std::string mismatch_string;
    int c = 0;
    if (offset >= 0) {
        size_t overlap = std::min(correct_seq.size() - std::abs(offset), mismatch_seq.size());
        for (size_t pos = 0; pos < overlap; ++pos) {
            if (mismatch_seq[pos] != correct_seq[pos + offset]) {
                if (c > 0) {
                    if (c == 1) mismatch_string += mismatch_seq[pos - 1];
                    else mismatch_string += std::to_string(c);
                }
                mismatch_string += mismatch_seq[pos];
                c = 0;
            } else {
                c++;
            }
        }
    } else {
        size_t overlap = std::min(mismatch_seq.size() - std::abs(offset), correct_seq.size());
        for (size_t pos = 0; pos < overlap; ++pos) {
            if (mismatch_seq[pos - offset] != correct_seq[pos]) { // 注意此时 offset<0
                if (c > 0) {
                    if (c == 1) mismatch_string += mismatch_seq[pos - offset - 1];
                    else mismatch_string += std::to_string(c);
                }
                mismatch_string += mismatch_seq[pos - offset];
                c = 0;
            } else {
                c++;
            }
        }
    }
    return mismatch_string;
}

std::string RefFreeSeqProcess::get_seq(uint32_t read_id) {
    std::string seq;
    seq.reserve(seq_view[read_id].size());
    if (is_contain_N[read_id]) {
        for (size_t i = 0; i < seq_view[read_id].size(); ++i) {
            if (N_encode[read_id][i] == '1') seq += 'N';
            else seq += seq_view[read_id][i];
        }
    } else {
        for (auto ch:seq_view[read_id]) seq += ch;
    }
    return seq;
}

void RefFreeSeqProcess::graph_encode(Graph &graph) {
    using namespace lemon;
    DEFINE_PATH_MACRO;

    std::ofstream extra_read_txt(extra_read_path.c_str());
    std::ofstream read_txt(read_path.c_str());
    std::ofstream overlap_dir_txt(overlap_dir_path.c_str());
    std::ofstream strand_dir_txt(strand_dir_path.c_str());
    std::ofstream mismatch_txt(mismatch_path.c_str());
    std::ofstream tree_txt(tree_path.c_str());

    Graph::G::NodeMap<u8> node_strand(graph.g, 0);
    Graph::G::ArcMap<bool> filter_edge_map(graph.g);
    std::vector<u32> extra_read_id_list;
    std::vector<u32> tree_read_id_list;
    std::vector<u32> overlap_offset;

    kruskal(graph.g, graph.weight_map, filter_edge_map);
    auto filter_dir_g = filterArcs(graph.g, filter_edge_map);
    using F_G = decltype(filter_dir_g);
    F_G::NodeMap<int> component_id(filter_dir_g);
    auto component_number = connectedComponents(undirector(filter_dir_g), component_id);
    std::vector<F_G::Node> component;
    component.resize(component_number, INVALID);
    for (F_G::NodeIt node(filter_dir_g); node != INVALID; ++node) {
        if (countOutArcs(filter_dir_g, node) == 0 && countInArcs(filter_dir_g, node) == 0) {
            extra_read_id_list.push_back(filter_dir_g.id(node));
            continue;
        }
        if (component[component_id[node]] == INVALID) {
            component[component_id[node]] = node;
        }
    }

    for (const auto &root : component) {
        if (root != INVALID) {
            auto tree_string = get_seq(filter_dir_g.id(root));
            read_txt << tree_string;
            node_strand[root] = 0;
            std::deque<std::pair<F_G::Node, F_G::Node>> q;
            q.emplace_back(root, root);
            F_G::Node last_p_node = root;

            while (!q.empty()) {
                // auto[n_node, p_node] = q.front();
                auto n_node = q.front().first;
                auto p_node = q.front().second;
                q.pop_front();

                u32 n_node_id = filter_dir_g.id(n_node);
                tree_read_id_list.push_back(n_node_id);
                std::string n_seq = get_seq(n_node_id);
                if (node_strand[n_node] == 1) reverse_complement(n_seq);
                strand_dir_txt << std::to_string(node_strand[n_node]);

                int child_num = 0;
                for (F_G::OutArcIt out_e(filter_dir_g, n_node); out_e != INVALID; ++out_e) {
                    auto child_node = filter_dir_g.target(out_e);
                    if (child_node != p_node) {
                        u32 child_node_id = filter_dir_g.id(child_node);

                        u8 overlap_dir = 0;
                        node_strand[child_node] = graph.strand_map[out_e].end_strand_id;
                        int off = graph.weight_map[out_e].off;

                        if (graph.strand_map[out_e].start_strand_id != node_strand[n_node]) {
                            node_strand[child_node] = !graph.strand_map[out_e].end_strand_id;
                            overlap_dir = 1;
                            off = -graph.weight_map[out_e].off;
                        }

                        std::string child_seq = get_seq(child_node_id);

                        if (node_strand[child_node] == 1) reverse_complement(child_seq);

                        overlap_dir_txt << std::to_string(overlap_dir);
                        overlap_offset.push_back(graph.weight_map[out_e].off);

                        if ((overlap_dir == 0 && std::abs(off) >= n_seq.size()) ||
                            (overlap_dir == 1 && std::abs(off) >= child_seq.size())) {
                            // 虽然前面记录的是统一序列方向后的数据,但是这里依然用原先的比对方向计算.
                            std::string n_seq_ = n_seq, child_seq_ = child_seq;
                            reverse_complement(n_seq_); reverse_complement(child_seq_);
                            read_txt << compute_overlap(child_seq_, n_seq_, !overlap_dir, graph.weight_map[out_e].off);
                            mismatch_txt << compute_mismatch_string(n_seq_, child_seq_, -off) << "\n";
                        } else {
                            // 正常执行统一序列方向后的比对
                            read_txt << compute_overlap(child_seq, n_seq, overlap_dir, graph.weight_map[out_e].off);
                            mismatch_txt << compute_mismatch_string(n_seq, child_seq, off) << "\n";
                        }

                        q.emplace_back(child_node, n_node);
                        child_num++; // child_num 应该在这里更新
                    }
                }

                for (F_G::InArcIt in_e(filter_dir_g, n_node); in_e != INVALID; ++in_e) {
                    auto child_node = filter_dir_g.source(in_e);
                    if (child_node != p_node) {
                        u32 child_node_id = filter_dir_g.id(child_node);

                        u8 overlap_dir = 1;
                        node_strand[child_node] = graph.strand_map[in_e].start_strand_id;
                        int off = -graph.weight_map[in_e].off;
                        if (graph.strand_map[in_e].end_strand_id != node_strand[n_node]) {
                            overlap_dir = 0;
                            node_strand[child_node] = !graph.strand_map[in_e].start_strand_id;
                            off = graph.weight_map[in_e].off;
                        }

                        std::string child_seq = get_seq(child_node_id);

                        if (node_strand[child_node] == 1) reverse_complement(child_seq);

                        overlap_dir_txt << std::to_string(overlap_dir);
                        overlap_offset.push_back(graph.weight_map[in_e].off);

                        if ((overlap_dir == 0 && std::abs(off) >= n_seq.size()) ||
                            (overlap_dir == 1 && std::abs(off) >= child_seq.size())) {
                            // 虽然前面记录的是统一序列方向后的数据,但是这里依然用原先的比对方向计算.
                            std::string n_seq_ = n_seq, child_seq_ = child_seq;
                            reverse_complement(n_seq_); reverse_complement(child_seq_);
                            read_txt << compute_overlap(child_seq_, n_seq_, !overlap_dir, graph.weight_map[in_e].off);
                            mismatch_txt << compute_mismatch_string(n_seq_, child_seq_, -off) << "\n";
                        } else {
                            // 正常执行统一序列方向后的比对
                            read_txt << compute_overlap(child_seq, n_seq, overlap_dir, graph.weight_map[in_e].off);
                            mismatch_txt << compute_mismatch_string(n_seq, child_seq, off) << "\n";
                        }

                        q.emplace_back(child_node, n_node);
                        child_num++; // child_num 应该在这里更新
                    }
                }

                // 根据父节点的改变,确定加上区间分割符"1"
                if (last_p_node != p_node) {
                    tree_txt << "1";
                    last_p_node = p_node;
                }

                if (n_node != p_node) { // 树结构编码字符串中,不记录根节点的信息
                    if (child_num > 0) {
                        tree_txt << "0";
                    } else { // 叶子节点
                        tree_txt << "2";
                    }
                }
            }
            tree_txt << "1"; // 最后加1分开不同的树
        }
    }

    {
        std::vector<char> overlap_offset_compress;
        Integer_compress(overlap_offset, overlap_offset_compress);
        write_vector_to_binary_file(overlap_offset_compress, overlap_offset_path.c_str());
    }

    std::sort(extra_read_id_list.begin(), extra_read_id_list.end());
    for (u32 extra_read_id : extra_read_id_list) {
        extra_read_txt << get_seq(extra_read_id);
    }

    if (m_param->m_seq_compress_mode == SEQ_REORDER) {
        uint32_t order_id = 0;
        order.resize(seq_view.size());
        for (u32 tree_read_id : tree_read_id_list) {
            order[order_id++] = tree_read_id; // 提供新的顺序ID
        }
        for (u32 extra_read_id : extra_read_id_list) {
            order[order_id++] = extra_read_id;
        }
        if (m_param->m_arcType == ARCTYPE_PE) {
            std::ofstream PE_file(pe_id_path.c_str());
            for (u32 tree_read_id : tree_read_id_list) PE_file << std::to_string(tree_read_id % 2);
            for (u32 extra_read_id : extra_read_id_list) PE_file << std::to_string(extra_read_id % 2);
        }
    } else {
        // 只需要保存tree_read_id_list. extra_read按照原始的顺序输出,等tree_read解压并填充结束后再顺序填充extra_read即可.
        std::vector<char> tree_read_id_list_compress;
        Integer_compress(tree_read_id_list, tree_read_id_list_compress);
        write_vector_to_binary_file(tree_read_id_list_compress, tree_read_id_path.c_str());
    }

}

int RefFreeSeqProcess::encode_output(char *output) {
    DEFINE_PATH_MACRO;
    if (!boost::filesystem::exists(parallel_command_path)) {
        throw std::runtime_error("PARALLEL_COMMAND_PATH isn't a valid path");
    }
    std::string bsc_compress_command = std::string(parallel_command_path) + " ";
    if (m_param->m_seq_compress_mode == SEQ_PRESERVE_ORDER) {
        bsc_compress_command += " \"" + bsc_compress(tree_read_id_path.c_str(), tree_read_id_comp_path.c_str()) + "\" ";
    } else if (m_param->m_arcType == ARCTYPE_PE) { // reorder && pe_file
        bsc_compress_command += " \"" + bsc_compress(pe_id_path.c_str(), pe_id_comp_path.c_str()) + "\" ";
    }
    bsc_compress_command += " \"" + bsc_compress(read_path.c_str(), read_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(overlap_dir_path.c_str(), overlap_dir_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(strand_dir_path.c_str(), strand_dir_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(mismatch_path.c_str(), mismatch_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(tree_path.c_str(), tree_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(extra_read_path.c_str(), extra_read_comp_path.c_str()) + "\" ";
    bsc_compress_command += " \"" + bsc_compress(overlap_offset_path.c_str(), overlap_offset_comp_path.c_str()) + "\" ";
    bsc_compress_command += " > /dev/null";  // suppress output
    int status = std::system(bsc_compress_command.c_str());
    if (status != 0) throw std::runtime_error("Error occurred during bsc compress.");

    // 编码输出,先处理顺序这一部分,然后按照 read overlap_dir strand_dir mismatch tree_path extra_read overlap_off 的顺序保存
    int compress_size = 0;
    auto read_file_and_output = [&](const std::string &filename) -> size_t {
        std::vector<uint8_t> data;
        read_vector_from_binary_file(data, filename);
        u32 size = data.size();
        std::memcpy(output, &size, sizeof(u32));    // write stream size
        output += sizeof(u32);
        std::memcpy(output, data.data(), size); // write stream
        output += size;
        return size + sizeof(u32);
    };
    if (m_param->m_seq_compress_mode == SEQ_PRESERVE_ORDER) {
        compress_size += read_file_and_output(tree_read_id_comp_path.c_str());
    } else if (m_param->m_arcType == ARCTYPE_PE) {
        compress_size += read_file_and_output(pe_id_comp_path.c_str());
    }
    compress_size += read_file_and_output(read_comp_path.c_str());
    compress_size += read_file_and_output(overlap_dir_comp_path.c_str());
    compress_size += read_file_and_output(strand_dir_comp_path.c_str());
    compress_size += read_file_and_output(mismatch_comp_path.c_str());
    compress_size += read_file_and_output(tree_comp_path.c_str());
    compress_size += read_file_and_output(extra_read_comp_path.c_str());
    compress_size += read_file_and_output(overlap_offset_comp_path.c_str());

    return compress_size;
}

int RefFreeSeqProcess::internal_compress(char *output) {
    Graph graph(seq_view.size());
    build_minimizer_index_and_alignment(graph);
    graph_encode(graph);
    auto compress_size = encode_output(output);
    return compress_size;
}

void RefFreeSeqProcess::clean_working_directory() {
    // for (const auto &entry : boost::filesystem::directory_iterator(working_path))
    //     boost::filesystem::remove_all(entry.path());
}

void RefFreeSeqProcess::clear_all() {
    max_seq_len = 0;
    is_long_read_compress = false;
    std::vector<nonstd::string_view>().swap(seq_view);
    std::vector<std::string>().swap(N_reconstruct);
    std::vector<std::string>().swap(N_encode);
    std::vector<bool>().swap(is_contain_N);
    std::vector<uint32_t>().swap(order);
    std::vector<uint8_t>().swap(pe_id);
}

int RefFreeSeqProcess::calcMd5(char *outptr) {
    if(m_paramPtr->m_seq_compress_mode == SEQ_PRESERVE_ORDER)
    {
        if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
        {
            calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &m_hashval);
            memcpy(outptr, &m_hashval, 16);
        }
        else if(m_paramPtr->m_checktype == DATACHECK_MD5)
        {
            MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), m_szmd5);
            memcpy(outptr, m_szmd5, 16);
        }
        return 16;
    }
    return 0;
}


int RefFreeSeqProcess::long_read_compress_implement(char *outptr) {
    return 0;
}

int RefFreeSeqProcess::long_read_compress(char *outptr) {
    int compress_size = long_read_compress_implement(outptr);
    if (compress_size == 0) {
        throw std::runtime_error("RefFree long read compress is unsupport now");
    }
    clean_working_directory();
    return compress_size;
}

int RefFreeSeqProcess::do_compress_seq(char *outptr) {
    clear_all();
    auto *seq_memory = m_seqElPtr->m_rawBuf.getdata();
    size_t seq_index = 0;

    seq_view.reserve(m_seqLenElPtr->m_LenArry.size());
    is_contain_N.reserve(m_seqLenElPtr->m_LenArry.size());

    for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); i++) 
    {
        seq_view.emplace_back(seq_memory + seq_index, m_seqLenElPtr->m_LenArry[i]);
        is_contain_N.push_back(seq_view.back().find('N') != nonstd::string_view::npos);
        max_seq_len = std::max(max_seq_len, m_seqLenElPtr->m_LenArry[i]);
        if (!is_long_read_compress && (m_seqLenElPtr->m_LenArry[i] > MAX_SHORT_READ_LEN)) 
        {
            is_long_read_compress = true;
        }
        seq_index += m_seqLenElPtr->m_LenArry[i];
    }
    if (is_long_read_compress) {
        return long_read_compress(outptr);
    }
    N_encode.resize(is_contain_N.size());
    N_reconstruct.resize(is_contain_N.size());
    for (size_t i = 0; i < is_contain_N.size(); ++i) {
        if (is_contain_N[i]) {
            Count count;
            for (auto ch:seq_view[i]) count.add(ch);
            auto most_frequent_base = count.get();
            for (auto ch:seq_view[i]) {
                if (ch == 'N') {
                    N_reconstruct[i] += most_frequent_base;
                    N_encode[i] += '1';
                } else {
                    N_reconstruct[i] += ch;
                    N_encode[i] += '0';
                }
            }
            seq_view[i] = N_reconstruct[i]; // seq_view 每条序列都是不含N的
        }
    }
    int compress_size = internal_compress(outptr);
    clean_working_directory();
    if (compress_size == 0) {
        arc_stdout(level::ERROR, "RefFreeSeqProcess compress error");
    }
    return compress_size;
}

int RefFreeSeqProcess::compress(char *outptr) {
    auto timer1 = m_profile->getTimer();
    timer1->start("SeqMd5Time");
    int szmd5 = calcMd5(outptr);
    outptr += szmd5;
    timer1->stop();

    timer1->start("SeqTime");
    int szseq = do_compress_seq(outptr);
    timer1->stop();
    return szmd5 + szseq;
}

bool RefFreeSeqProcess::checkData() {
    if(m_paramPtr->m_seq_compress_mode == SEQ_PRESERVE_ORDER)
    {
        if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
        {
            XXH128_hash_t res;
            bool ret = calcDataHash(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), &res);
            if(!ret || XXH128_isEqual(res, m_hashval)==0)
            {
                std::ofstream debug_seq((working_path/"debug_seq.txt").c_str());
                auto start = m_seqElPtr->m_rawBuf.getdata();
                for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); ++i) {
                    debug_seq << std::string(start, m_seqLenElPtr->m_LenArry[i]) << "\n" ;
                    start += m_seqLenElPtr->m_LenArry[i];
                }
                arc_stdout(level::ERROR, "seq md5 check fail");
            }
        }
        else if(m_paramPtr->m_checktype == DATACHECK_MD5)
        {
            unsigned char szmd5[16] = {0};
            MDString(m_seqElPtr->m_rawBuf.getdata(), m_seqElPtr->m_rawBuf.size(), szmd5);

            if(memcmp(m_szmd5, szmd5, 16) != 0)
            {
                std::ofstream debug_seq((working_path/"debug_seq.txt").c_str());
                auto start = m_seqElPtr->m_rawBuf.getdata();
                for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); ++i) {
                    debug_seq << std::string(start, m_seqLenElPtr->m_LenArry[i]) << "\n" ;
                    start += m_seqLenElPtr->m_LenArry[i];
                }
                arc_stdout(level::ERROR, "seq md5 check fail");
            }
        }
    }

    return true;
}

int RefFreeSeqProcess::decompress(int inlen, char *inptr) {
    if(m_paramPtr->m_seq_compress_mode == SEQ_PRESERVE_ORDER)
    {
        if(m_paramPtr->m_checktype == DATACHECK_XXHASH)
        {
            memcpy(&m_hashval, inptr, 16);
        }
        else if(m_paramPtr->m_checktype == DATACHECK_MD5)
        {
            memcpy(m_szmd5, inptr, 16);
        }
        inptr += 16;
        inlen -= 16;
    }

    clear_all();
    DEFINE_PATH_MACRO;

    std::vector<u32> tree_read_id_list;
    std::vector<u32> overlap_offset;

    auto read_buffer_to_file = [&](const boost::filesystem::path &bsc_filepath,
                                   const boost::filesystem::path &bsc_decompress_filepath) -> size_t {
        uint32_t buffer_size;
        std::memcpy(&buffer_size, inptr, sizeof(uint32_t));
        inptr += sizeof(uint32_t);
        {
            std::ofstream bsc_file(bsc_filepath.c_str(), std::ios::binary);
            bsc_file.write(inptr, buffer_size);
            bsc_file.flush();
            inptr += buffer_size;
        } // auto close bsc_file after leaving scope
        bsc_decompress(bsc_filepath.c_str(), bsc_decompress_filepath.c_str());
        return sizeof(uint32_t) + buffer_size;
    };
    size_t total_read_size = 0;

    if (m_param->m_seq_compress_mode == SEQ_PRESERVE_ORDER) { // 保留顺序
        total_read_size += read_buffer_to_file(tree_read_id_comp_path, tree_read_id_path);
        std::vector<char> tree_read_id_list_compress;
        read_vector_from_binary_file(tree_read_id_list_compress, tree_read_id_path.c_str());
        Integer_decompress(tree_read_id_list_compress, tree_read_id_list);
    }
    if (m_param->m_seq_compress_mode == SEQ_REORDER && m_param->m_arcType == ARCTYPE_PE) { // 不保留顺序但为PE文件
        total_read_size += read_buffer_to_file(pe_id_comp_path, pe_id_path);
        std::ifstream pe_id_file(pe_id_path.c_str());
        char c;
        pe_id.reserve(m_seqLenElPtr->m_LenArry.size());
        while (pe_id_file >> c) {
            pe_id.push_back(c - '0'); // 0 or 1
        }
    }

    total_read_size += read_buffer_to_file(read_comp_path, read_path);
    total_read_size += read_buffer_to_file(overlap_dir_comp_path, overlap_dir_path);
    total_read_size += read_buffer_to_file(strand_dir_comp_path, strand_dir_path);
    total_read_size += read_buffer_to_file(mismatch_comp_path, mismatch_path);
    total_read_size += read_buffer_to_file(tree_comp_path, tree_path);
    total_read_size += read_buffer_to_file(extra_read_comp_path, extra_read_path);
    total_read_size += read_buffer_to_file(overlap_offset_comp_path, overlap_offset_path);
    {
        std::vector<char> overlap_offset_compress;
        read_vector_from_binary_file(overlap_offset_compress, overlap_offset_path.c_str());
        Integer_decompress(overlap_offset_compress, overlap_offset);
    }

    if (total_read_size != inlen) {
        arc_stdout(level::ERROR, "Read Length not equal to in length");
    }

    std::ifstream extra_read_txt(extra_read_path.c_str());
    std::ifstream read_txt(read_path.c_str());
    std::ifstream overlap_dir_txt(overlap_dir_path.c_str());
    std::ifstream strand_dir_txt(strand_dir_path.c_str());
    std::ifstream mismatch_txt(mismatch_path.c_str());
    std::ifstream tree_txt(tree_path.c_str());

    std::vector<int> seq_dir(m_seqLenElPtr->m_LenArry.size(),
                             -1); /// seq_dir : -1 代表未解析的序列; 0 代表储存forward序列; 1代表储存rc序列(reverse complement)
    std::vector<char *> sequences(m_seqLenElPtr->m_LenArry.size() + 1); /// 直接对要解压的数组进行操作

    size_t capacity = m_seqElPtr->m_rawBuf.capacity();
    size_t decompress_size = 0;
    sequences[0] = (char *) m_seqElPtr->m_rawBuf.getdata();
    for (size_t i = 0; i < m_seqLenElPtr->m_LenArry.size(); ++i) {
        uint32_t seq_len = m_seqLenElPtr->m_LenArry[i];
        decompress_size += seq_len;
        sequences[i + 1] = sequences[i] + seq_len;
    }
    if (decompress_size > capacity) {
        arc_stdout(level::ERROR, "Decompress size larger than capacity");
        return 0;
    }
    m_seqElPtr->m_rawBuf.updatesize(decompress_size);

    uint32_t read_index = 0;
    uint32_t overlap_offset_index = 0;

    auto get_read_id = [&]() -> uint32_t {
        if (m_param->m_seq_compress_mode == SEQ_PRESERVE_ORDER) {
            return tree_read_id_list[read_index++];
        }
        return read_index++;
    };
    size_t size = boost::filesystem::file_size(read_path);
    size_t process_cnt = 0;
    while (process_cnt < size) {
        auto read_id = get_read_id();
        auto read_len = m_seqLenElPtr->m_LenArry[read_id];
        for (size_t i = 0; i < read_len; ++i) {
            char ch;
            read_txt >> noskipws >> ch;
            process_cnt++;
            sequences[read_id][i] = ch;
        }
        char strand_dir;
        strand_dir_txt >> strand_dir;
        seq_dir[read_id] = strand_dir - '0';
        std::queue<uint32_t> Q;
        Q.push(read_id);
        while (!Q.empty()) {
            auto parent_id = Q.front();
            auto parent_len = m_seqLenElPtr->m_LenArry[parent_id];
            auto parent_str = sequences[parent_id];
            Q.pop();

            char tree_op;
            while (true) {
                tree_txt >> tree_op;
                if (tree_op == '1') break; /// 分隔符
                auto child_id = get_read_id();
                auto child_str = sequences[child_id];
                auto child_len = m_seqLenElPtr->m_LenArry[child_id];
                strand_dir_txt >> strand_dir;
                seq_dir[child_id] = strand_dir - '0';
                char overlap_dir;
                overlap_dir_txt >> overlap_dir;
                std::string mismatch_str;
                std::unordered_map<int, char> mismatch_map;
                std::getline(mismatch_txt, mismatch_str);
                int pos = 0;
                std::string num = "0";
                for (auto ch: mismatch_str) {
                    if ('0' <= ch && ch <= '9') num += ch;
                    else {
                        pos += std::stoi(num);
                        mismatch_map.emplace(pos, ch);
                        num = "0";
                        pos++;
                    }
                }
                int off = overlap_offset[overlap_offset_index++];

                if ((overlap_dir == '0' && off >= parent_len) ||
                    (overlap_dir == '1' && off >= child_len) ){
                    // 反过来解压序列
                    std::string parent_str_ (parent_str, parent_len);
                    reverse_complement(parent_str_);
                    if (overlap_dir == '1') { // 反过来是后缀前缀重叠
                        int i = 0;
                        size_t overlap_range = std::min(parent_len - off, child_len);
                        for (; i < overlap_range; ++i) {
                            auto itr = mismatch_map.find(i);
                            if (itr != mismatch_map.end()) {
                                child_str[i] = itr->second;
                            } else {
                                child_str[i] = parent_str_[i + off];
                            }
                        }
                        int leave_size = off + (int) child_len - (int) parent_len;
                        if (leave_size > 0) {
                            while (leave_size--) {
                                char ch;
                                read_txt >> noskipws >> ch;
                                process_cnt++;
                                child_str[i++] = ch;
                            }
                        }
                    } else {  // 反过来是前缀后缀重叠
                        for (int i = 0; i < off; ++i) {
                            char ch;
                            read_txt >> noskipws >> ch;
                            process_cnt++;
                            child_str[i] = ch;
                        }
                        int i = 0;
                        size_t overlap_range = std::min(child_len - off, parent_len);
                        for (; i < overlap_range; ++i) {
                            auto itr = mismatch_map.find(i);
                            if (itr != mismatch_map.end()) {
                                child_str[i + off] = itr->second;
                            } else {
                                child_str[i + off] = parent_str_[i];
                            }
                        }
                        int leave_size = (int) child_len - (off + (int) parent_len);
                        if (leave_size > 0) {
                            while (leave_size--) {
                                char ch;
                                read_txt >> noskipws >> ch;
                                process_cnt++;
                                child_str[off + i] = ch;
                                i++;
                            }
                        }
                    }
                    reverse_complement(child_str, child_len); // 最后再将child_str方向反过来
                } else {
                    if (overlap_dir == '0') { // 后缀-前缀重叠
                        int i = 0;
                        size_t overlap_range = std::min(parent_len - off, child_len);
                        for (; i < overlap_range; ++i) {
                            auto itr = mismatch_map.find(i);
                            if (itr != mismatch_map.end()) {
                                child_str[i] = itr->second;
                            } else {
                                child_str[i] = parent_str[i + off];
                            }
                        }
                        int leave_size = off + (int) child_len - (int) parent_len;
                        if (leave_size > 0) {
                            while (leave_size--) {
                                char ch;
                                read_txt >> noskipws >> ch;
                                process_cnt++;
                                child_str[i++] = ch;
                            }
                        }
                    } else { // 前缀-后缀重叠
                        for (int i = 0; i < off; ++i) {
                            char ch;
                            read_txt >> noskipws >> ch;
                            process_cnt++;
                            child_str[i] = ch;
                        }
                        int i = 0;
                        size_t overlap_range = std::min(child_len - off, parent_len);
                        for (; i < overlap_range; ++i) {
                            auto itr = mismatch_map.find(i);
                            if (itr != mismatch_map.end()) {
                                child_str[i + off] = itr->second;
                            } else {
                                child_str[i + off] = parent_str[i];
                            }
                        }
                        int leave_size = (int) child_len - (off + (int) parent_len);
                        if (leave_size > 0) {
                            while (leave_size--) {
                                char ch;
                                read_txt >> noskipws >> ch;
                                process_cnt++;
                                child_str[off + i] = ch;
                                i++;
                            }
                        }
                    }
                }
                if (tree_op == '0') {
                    Q.push(child_id);
                } /// 如果 tree_op == '2', 那么已经是叶子节点,不需要继续处理...
            }
        }
    }

    for (size_t seq_id = 0; seq_id < m_seqLenElPtr->m_LenArry.size(); ++seq_id) {
        if (seq_dir[seq_id] == -1) {   // 还原extra_read
            auto seq = sequences[seq_id];
            int seq_len = m_seqLenElPtr->m_LenArry[seq_id];
            for (int i = 0; i < seq_len; ++i) {
                char ch;
                extra_read_txt >> noskipws >> ch;
                seq[i] = ch;
            }
        } else if (seq_dir[seq_id] == 1) { // 从rc方向改成forward方向
            auto seq = sequences[seq_id];
            int seq_len = m_seqLenElPtr->m_LenArry[seq_id];
            reverse_complement(seq, seq_len);
        }
    }

    checkData();
    clean_working_directory();
    return decompress_size;
}