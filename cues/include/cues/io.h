/**
 * Lightweight wrappers for HtsLib-based FASTA and BAM IO
 */
#ifndef CUES_INCLUDE_CUES_IO_H_
#define CUES_INCLUDE_CUES_IO_H_

#include <algorithm>
#include <cctype>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/kseq.h>
#include <htslib/sam.h>
#include <iostream>
#include <filesystem>
#include <format>
#include <optional>
#include <ranges>
#include <string_view>
#include <tuple>
#include <vector>
#include <utility>
#include <unordered_map>
#include <zlib.h>

namespace cues::io {

    /////////////////
    ///  FASTA IO ///
    /////////////////
    class Fasta {
    public:
        explicit Fasta(const std::string &fasta_file_name) {
            if (fasta_file_name.empty()) {
                throw std::invalid_argument("[FASTA IO] Reference fasta file was not provided (empty).");
            }
            fai_ = fai_load(fasta_file_name.c_str());
            if (!fai_) {
                throw std::runtime_error("[FASTA IO] Failed to open FASTA file: " + fasta_file_name);
            }
        }
        Fasta(const Fasta&) = delete;
        Fasta& operator=(const Fasta&) = delete;
        Fasta(Fasta&& other) noexcept : fai_(std::exchange(other.fai_, nullptr)) {}
        Fasta& operator=(Fasta&& other) noexcept {
            if (this != &other) {
                if (fai_) fai_destroy(fai_);
                fai_ = std::exchange(other.fai_, nullptr);
            }
            return *this;
        }
        ~Fasta() {
            if (fai_) {
                fai_destroy(fai_);
            }
        }

        inline std::string fetch(const std::string& chr_name,
                                 std::optional<hts_pos_t> start = std::nullopt,
                                 std::optional<hts_pos_t> end = std::nullopt) {
            int chr_len = faidx_seq_len(fai_, chr_name.c_str());
            if (chr_len < 0) throw std::runtime_error("[FASTA IO] Missing chromosome: " + chr_name);
            hts_pos_t ref_start = start.value_or(0);
            hts_pos_t ref_end = end.value_or(static_cast<hts_pos_t>(chr_len));
            if (ref_start >= ref_end || ref_end < 1 || ref_end > chr_len) {
                throw std::runtime_error(std::format("[FASTA IO] Invalid reference range {}: {}-{}", chr_name, ref_start, ref_end));
            }
            hts_pos_t fetched_seq_len;
            char* seq = faidx_fetch_seq64(fai_, chr_name.c_str(), ref_start, ref_end - 1, &fetched_seq_len);
            std::string ref_seq(seq, fetched_seq_len);
            free(seq);
            return ref_seq;
        }

        inline std::vector<std::string> get_chr_names() const {
            const int n_seqs = faidx_nseq(fai_);
            std::vector<std::string> chr_names(n_seqs);
            for (int i = 0; i < n_seqs; i++) {
                const char* chr_name = faidx_iseq(fai_, i);
                chr_names[i] = chr_name;
            }
            return chr_names;
        }

        inline int64_t get_chr_len(const std::string& chr_name) const {
            return faidx_seq_len64(fai_, chr_name.c_str());
        }

    private:
        faidx_t* fai_;
    };

    /////////////////
    ///  FASTQ IO ///
    /////////////////
    KSEQ_INIT(gzFile, gzread)
    struct FastqReadDelete {
        void operator()(kseq_t* r) { kseq_destroy(r); }
    };
    using fastq_read_ptr_t = std::unique_ptr<kseq_t, FastqReadDelete>;
    struct FastqRead {
        std::string name;
        std::string seq;
    };
    using read_buffer_t = std::vector<FastqRead>;

    class Fastq {
    public:
        explicit Fastq(const std::string &fastq_file_name, const std::optional<int> batch_size = std::nullopt,
                       const std::optional<int> max_reads = std::nullopt):
                       batch_size_(batch_size), n_reads_processed_(0), max_reads_to_process_(max_reads)  {
            if (fastq_file_name.empty()) {
                throw std::invalid_argument("[FASTQ IO] FASTQ file was not provided (empty).");
            }
            fp_ = gzopen(fastq_file_name.c_str(), "r");
            if (!fp_) {
                throw std::runtime_error("[FASTQ IO] Failed to open FASTQ file: " + fastq_file_name);
            }
            seq_ = fastq_read_ptr_t(kseq_init(fp_));
            gzrewind(fp_);
            kseq_rewind(seq_);
            if (batch_size_.has_value()) {
                buffer_.resize(batch_size_.value());
            }
        }

        Fastq(const Fastq&) = delete;
        Fastq& operator=(const Fastq&) = delete;

        ~Fastq() {
            if (fp_) {
                gzclose(fp_);
                fp_ = NULL;
            }
        }

        [[nodiscard]] inline size_t get_n_processed() const { return n_reads_processed_; }
        [[nodiscard]] inline bool limit_reached() const {
            return max_reads_to_process_.has_value() && n_reads_processed_ >= max_reads_to_process_;
        }

        inline read_buffer_t& read_next_batch() {
            if (!batch_size_.has_value()) {
                throw std::runtime_error("[FASTQ IO] Batched-mode reading was not initialized.");
            }
            const auto processed_at_start = n_reads_processed_;
            for (auto i = 0; i < batch_size_.value(); i++) {
                if (limit_reached()) break; // reached the specified number of records
                int ret = kseq_read(seq_.get());
                if (ret < 0) {
                    if (ret == -1) break;
                    throw std::runtime_error("[FASTQ IO] Error reading FASTQ record");
                }
                buffer_[i].name = seq_->name.s;
                buffer_[i].seq = seq_->seq.s;
                n_reads_processed_++;
            }
            const auto n_processed = n_reads_processed_ - processed_at_start;
            if (n_processed < batch_size_.value()) {
                buffer_.resize(n_processed);
            }
            return buffer_;
        }

    private:
        gzFile fp_;
        fastq_read_ptr_t seq_;
        read_buffer_t buffer_;
        size_t n_reads_processed_;
        std::optional<size_t> batch_size_;
        std::optional<size_t> max_reads_to_process_;
    };


    /////////////////
    ///   BAM IO  ///
    /////////////////

    struct AlnDelete {
        void operator()(bam1_t *r) { bam_destroy1(r); }
    };
    using aln_ptr_t = std::unique_ptr<bam1_t, AlnDelete>;

    struct BAMHeaderDelete {
        void operator()(bam_hdr_t *h) { bam_hdr_destroy(h); }
    };
    struct SAMHeaderDelete {
        void operator()(sam_hdr_t *h) { sam_hdr_destroy(h); }
    };
    using bam_header_t = std::unique_ptr<bam_hdr_t, BAMHeaderDelete>;

    struct SplitAln {
        std::string read_name;
        std::string chr_name;
        int64_t ref_start;
        int64_t ref_end;
        int32_t query_start;
        int32_t query_end;
        bool is_reverse;
        uint16_t mapq;

        SplitAln(const std::string read_name, const std::string chr_name, const int64_t ref_start, const int64_t ref_end,
                 const int32_t query_start, const int32_t query_end, const bool is_reverse, const uint16_t mapq):
                 read_name(read_name), chr_name(chr_name), ref_start(ref_start), ref_end(ref_end),
                 query_start(query_start), query_end(query_end), is_reverse(is_reverse), mapq(mapq) {}
        SplitAln(const std::string read_name, const std::string chr_name, const int64_t ref_start, const bool is_reverse,
                 const uint16_t mapq, const std::string& cigar_str):
                 read_name(read_name), chr_name(chr_name), ref_start(ref_start), is_reverse(is_reverse), mapq(mapq) {
            uint32_t* cigar_buffer = nullptr;
            size_t mem = 0;
            char* end = nullptr;
            ssize_t n_cigar = sam_parse_cigar(cigar_str.c_str(), &end, &cigar_buffer, &mem);
            if (n_cigar < 0 || cigar_buffer == nullptr || end == cigar_str.c_str()) {
                free(cigar_buffer);
                throw std::runtime_error("[BAM IO] Failed to parse SA CIGAR: " + cigar_str);
            }
            std::vector<uint32_t> cigar(cigar_buffer, cigar_buffer + n_cigar);
            free(cigar_buffer);
            ref_end = ref_start + bam_cigar2rlen(cigar.size(), cigar.data()) - 1;
            query_start = 0;
            int query_aln_len = bam_cigar2qlen(cigar.size(), cigar.data());
            const int cigar_start_idx = is_reverse ? cigar.size() - 1 : 0;
            for (auto idx : std::vector{0, cigar.size()-1}) {
                if (bam_cigar_op(cigar[idx]) == BAM_CSOFT_CLIP ||
                    bam_cigar_op(cigar[idx]) == BAM_CHARD_CLIP) {
                    if (idx == cigar_start_idx) {
                        query_start += bam_cigar_oplen(cigar[idx]);
                    }
                    query_aln_len -= bam_cigar_oplen(cigar[idx]);
                }
            }
            query_end = query_aln_len + query_start;
        }
    };

    class BamAln {
    public:
        explicit BamAln() {
            aln = aln_ptr_t(bam_init1());
        }

        explicit BamAln(std::string chr_name): chr_name(chr_name) {
            aln = aln_ptr_t(bam_init1());
        }

        std::vector<std::tuple<hts_pos_t, hts_pos_t, hts_pos_t, bool>> cins; // cigar insertions (if parsed)

        inline void set_chr_name(bam_header_t& header) {
            if (!(aln->core.flag & BAM_FUNMAP)) {
                chr_name = header->target_name[aln->core.tid];
            }
        }
        inline void clear() {
            seq.clear();
            cins.clear();
        }

        inline SplitAln make_split() const {
            return SplitAln(get_name(), chr_name, get_ref_start(), get_ref_end(), get_query_start(),
                           get_query_start() + get_query_aln_len(), is_reverse(), get_mapq());
        }

        inline const std::string& get_chr_name() const { return chr_name; }
        inline bam1_t* get_aln_internal() const { return aln.get(); }
        inline std::string get_name() const { return bam_get_qname(aln.get()); }
        inline bool is_reverse() const { return (aln->core.flag & BAM_FREVERSE) != 0; }
        inline bool is_split() const { return has_tag("SA"); }
        inline uint16_t get_flag() const { return aln->core.flag; }
        inline uint16_t get_mapq() const { return aln->core.qual; }
        inline int64_t get_ref_start() const { return aln->core.pos; }
        inline int64_t get_ref_end() const { return bam_endpos(aln.get()) - 1; }
        inline uint32_t* get_cigar() const { return bam_get_cigar(aln.get()); }
        inline uint32_t get_cigar_len() const { return aln->core.n_cigar; }
        inline std::tuple<hts_pos_t, hts_pos_t, hts_pos_t> parse_cigar_in_range(const hts_pos_t ref_min,
                                                                           const hts_pos_t ref_max) const {
            hts_pos_t cins_total = 0;
            hts_pos_t cdel_total = 0;
            hts_pos_t cdiff_total = 0;
            hts_pos_t ref_pos = get_ref_start();
            int read_pos = 0;
            const uint32_t *cigar = get_cigar();
            for (auto i = 0; i < get_cigar_len(); i++) {
                if (ref_pos >= ref_max) break;
                const int op = bam_cigar_op(cigar[i]);
                const int op_len = bam_cigar_oplen(cigar[i]);
                if (op == BAM_CINS) {
                    if (ref_pos >= ref_min && ref_pos < ref_max) cins_total += op_len;
                    read_pos += op_len;
                } else if (op == BAM_CDEL) {
                    cdel_total += std::max(hts_pos_t(0), std::min(ref_max, ref_pos+op_len) - std::max(ref_min, ref_pos));
                    ref_pos += op_len;
                } else if (op == BAM_CDIFF) {
                    if (ref_pos >= ref_min && ref_pos < ref_max) cdiff_total += op_len;
                    ref_pos += op_len; read_pos += op_len;
                } else if (op == BAM_CMATCH || op == BAM_CEQUAL) {
                    ref_pos += op_len; read_pos += op_len;
                } else if (op == BAM_CHARD_CLIP || op == BAM_CSOFT_CLIP) {
                    read_pos += op_len;
                }
            }
            return {cins_total, cdel_total, cdiff_total};
        }

        inline double get_discordance() const {
            int n_match = 0;
            int n_mismatch = 0;
            const uint32_t *cigar = get_cigar();
            for (auto i = 0; i < get_cigar_len(); i++) {
                const int op = bam_cigar_op(cigar[i]);
                const int op_len = bam_cigar_oplen(cigar[i]);
                if (op == BAM_CDIFF || op == BAM_CDEL || op == BAM_CINS) n_mismatch += 1;
                if (op == BAM_CMATCH || op == BAM_CEQUAL) n_match += op_len;
            }
            return double(n_mismatch) / double(n_mismatch + n_match);
        }

        inline std::string& get_seq() {
            if (seq.empty()) set_seq();
            return seq;
        }
        inline void set_seq() {
            seq.clear();
            uint8_t *rseq = bam_get_seq(aln.get());
            for (int i = 0; i < aln->core.l_qseq; ++i) {
                seq += seq_nt16_str[bam_seqi(rseq, i)];
            }
        }
        inline int32_t get_query_start() const {
            const uint32_t *cigar = get_cigar();
            const int cigar_idx = is_reverse() ? get_cigar_len() - 1 : 0;
            if (cigar_idx >= 0 &&
               (bam_cigar_op(cigar[cigar_idx]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[cigar_idx]) == BAM_CHARD_CLIP)) {
                return bam_cigar_oplen(cigar[cigar_idx]);
            }
            return 0;
        }
        inline bool is_clipped() const {
            return get_clipping("left") > 0 or  get_clipping("right") > 0;
        }

        inline bool is_hard_clipped() const {
            const uint32_t *cigar = get_cigar();
            return bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP || bam_cigar_op(cigar[get_cigar_len() - 1]) == BAM_CHARD_CLIP;
        }

        inline int get_clipping(const std::string side) const {
            const uint32_t *cigar = get_cigar();
            const int cigar_idx = side == "right" ? get_cigar_len() - 1 : 0;
            if (bam_cigar_op(cigar[cigar_idx]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[cigar_idx]) == BAM_CHARD_CLIP) {
                return bam_cigar_oplen(cigar[cigar_idx]);
            }
            return 0;
        }

        inline int32_t get_query_len() const { return bam_cigar2qlen(get_cigar_len(), get_cigar()); }
        inline int32_t get_query_aln_len() const {
            int32_t query_aln_len = get_query_len();
            const uint32_t *cigar = get_cigar();
            for (auto cigar_idx : std::vector{0, get_cigar_len()-1}) {
                if (bam_cigar_op(cigar[cigar_idx]) == BAM_CSOFT_CLIP ||
                    bam_cigar_op(cigar[cigar_idx]) == BAM_CHARD_CLIP) {
                    query_aln_len -= bam_cigar_oplen(cigar[cigar_idx]);
                }
            }
            return query_aln_len;
        }
        inline bool has_tag(const std::string &tag_name) const {
            return bam_aux_get(aln.get(), tag_name.c_str()) != nullptr;
        }
        inline uint8_t *get_tag_data(const std::string &tag_name) const {
            uint8_t *tag_data = bam_aux_get(aln.get(), tag_name.c_str());
            if (!tag_data) throw std::runtime_error("[BAM IO] Missing BAM tag: " + tag_name);
            return tag_data;
        }
        inline std::vector<SplitAln> get_supplementary_alns() const {
            std::vector<SplitAln> supp;
            if (!has_tag("SA")) return supp;
            auto split_str = [](std::string& s, char delim) {
                return s | std::views::split(delim) | std::views::transform([](auto&& sub) {
                    return std::string(std::ranges::begin(sub), std::ranges::end(sub));
                });
            };
            std::string sa{bam_aux2Z(get_tag_data("SA"))};
            if (!sa.empty() && sa.ends_with(';')) sa.pop_back();
            for (auto sa_entry : split_str(sa, ';')) {
                auto sa_entry_fields = split_str(sa_entry, ',');
                std::vector<std::string> sa_fields(sa_entry_fields.begin(), sa_entry_fields.end());
                supp.push_back(SplitAln(get_name(), sa_fields[0], std::stoi(sa_fields[1])-1, sa_fields[2] == "-",
                                       std::stoi(sa_fields[4]), sa_fields[3]));
            }
            return supp;
        }
        template <typename T>
        inline T get_tag(const std::string &tag_name) const;
        template<typename T>
        inline T get_tag_array(const std::string &tag_name, int array_idx) const;

    private:
        std::string chr_name;
        std::string seq;
        aln_ptr_t aln;
    };

    template<>
    inline double BamAln::get_tag<double>(const std::string &tag_name) const { return bam_aux2f(get_tag_data(tag_name));}
    template<>
    inline int BamAln::get_tag<int>(const std::string &tag_name) const { return bam_aux2i(get_tag_data(tag_name)); }
    template<>
    inline int BamAln::get_tag_array<int>(const std::string &tag_name, const int array_idx) const {
        return bam_auxB2i(get_tag_data(tag_name), array_idx);
    }
    inline void revcomp(std::string& seq, std::string name) {
        const std::unordered_map<char, char> kDNAComplement = {
            {'A', 'T'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'}};
        std::reverse(seq.begin(), seq.end());
        for (char& base : seq) {
            if (kDNAComplement.contains(base)) {
                base = kDNAComplement.at(std::toupper(base));
            }
        }
    }
    inline std::string make_aln_uid(const SplitAln& aln) {
        return aln.read_name + ":" + std::to_string(aln.query_start);
    }
    inline std::string make_aln_uid(const BamAln& aln) {
        return aln.get_name() + ":" + std::to_string(aln.get_query_start());
    }

    using aln_buffer_t = std::vector<BamAln>;
    class BamFile {
    public:
        explicit BamFile(const std::string &bam_file_name, const std::string &mode) {
            if (bam_file_name.empty()) {
                throw std::invalid_argument("[IO] BAM file was not provided (empty).");
            }
            fp_ = sam_open(bam_file_name.c_str(), mode.c_str());
            if(!fp_) {
                throw std::runtime_error("[IO] Failed to open BAM file: " + bam_file_name);
            }
            idx_ = sam_index_load(fp_, bam_file_name.c_str());
        }

        ~BamFile() {
            if (fp_) sam_close(fp_);
            fp_ = NULL;
        }

        void close() {
            if (fp_) sam_close(fp_);
            fp_ = NULL;
        }

    protected:
        samFile *fp_;
        hts_idx_t *idx_;
    };

    class BamWriter : public BamFile {
    public:
        BamWriter(const std::string &bam_file_name, const bam_header_t &header) :
                BamFile(bam_file_name, "wb"), header_(header) {
            if (sam_hdr_write(fp_, header_.get()) < 0) {
                throw std::runtime_error("[IO] Failed to write BAM header to output file.");
            }
        }

        void write(const aln_ptr_t& read) {
            if (sam_write1(fp_, header_.get(), read.get()) < 0) {
                throw std::runtime_error("[IO] Failed to write read to output BAM file.");
            }
        }

    private:
        const bam_header_t &header_;
    };

    class BamReader : public BamFile {
    public:
        explicit BamReader(const std::string &bam_file_name,
                           const std::optional<int> batch_size = std::nullopt,
                           const std::optional<int> max_reads = std::nullopt) :
                BamFile(bam_file_name, "r"), batch_size_(batch_size), n_reads_processed_(0),
                max_reads_to_process_(max_reads) {
            header_ = bam_header_t(sam_hdr_read(fp_));
            if (!header_) {
                throw std::runtime_error("[IO] Failed to read BAM header.");
            }
            if (batch_size_.has_value()) {
                init_buffer();
            }
        }

        [[nodiscard]] inline bam_header_t const &get_header() const { return header_; }
        [[nodiscard]] inline int get_chr_id(const std::string& chr_name) const {
            const auto tid = bam_name2id(header_.get(), chr_name.c_str());
            if (tid < 0) {
                throw std::runtime_error("[IO] Chromosome name not in the BAM header: " + chr_name);
            }
            return tid;
        }
        [[nodiscard]] inline size_t get_n_processed() const { return n_reads_processed_; }
        [[nodiscard]] inline bool limit_reached() const {
            return max_reads_to_process_.has_value() && n_reads_processed_ >= max_reads_to_process_;
        }

        aln_buffer_t read_region_overlap(const std::string& chr_name, const int64_t ref_start, const int64_t ref_end,
                                         const std::optional<int> max_aln_fetched = std::nullopt) {
            if (!idx_) throw std::runtime_error("BAM index is required for this operation but is missing/failed to load");
            hts_itr_t* iter = sam_itr_queryi(idx_, get_chr_id(chr_name), ref_start, ref_end);
            if (!iter) {
                throw std::runtime_error("[IO] Could not fetch region " + std::format("{}:{}-{}", chr_name, ref_start, ref_end));
            }
            aln_buffer_t alns_at_region;
            auto aln = BamAln(chr_name);
            while (sam_itr_next(fp_, iter, aln.get_aln_internal()) >= 0) {
                alns_at_region.push_back(std::move(aln));
                aln = BamAln(chr_name);
                if (max_aln_fetched.has_value() and alns_at_region.size() >= max_aln_fetched.value()) break;
            }
            hts_itr_destroy(iter);
            return alns_at_region;
        }

        aln_buffer_t &read_next_batch() {
            if (!batch_size_.has_value()) {
                throw std::runtime_error("[IO] Batched-mode reading was not initialized.");
            }
            const auto processed_at_start = n_reads_processed_;
            for (auto i = 0; i < batch_size_.value(); i++) {
                if (limit_reached()) break; // reached the specified number of records
                buffer_[i].clear();
                if (sam_read1(fp_, header_.get(), buffer_[i].get_aln_internal()) < 0) break;
                buffer_[i].set_chr_name(header_);
                n_reads_processed_++;
            }
            const auto n_processed = n_reads_processed_ - processed_at_start;
            if (n_processed < batch_size_.value()) {
                buffer_.resize(n_processed);
            }
            return buffer_;
        }

        std::optional<BamAln> read() {
            BamAln aln;
            const int ret = sam_read1(fp_, header_.get(), aln.get_aln_internal());
            if (ret < 0) return std::nullopt;
            aln.set_chr_name(header_);
            return aln;
        }

    private:
        void init_buffer() {
            buffer_.resize(batch_size_.value());
            for (int i = 0; i < batch_size_; i++) {
                buffer_[i] = BamAln();
            }
        }

        bam_header_t header_;
        aln_buffer_t buffer_;
        size_t n_reads_processed_;
        std::optional<size_t> batch_size_;
        std::optional<size_t> max_reads_to_process_;
    };

#define CUES_LOG(msg) do { \
    std::cout << msg << std::endl; \
} while (0)

}

#endif //CUES_INCLUDE_CUES_IO_H_
