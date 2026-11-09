/**
 * Lightweight wrapper around minimap2 for read alignment
 */
#ifndef CUES_INCLUDE_CUES_ALIGN_H_
#define CUES_INCLUDE_CUES_ALIGN_H_

#include "cues/io.h"
#include <format>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include <iostream>
#include <filesystem>
#include <minimap.h>
#include <vector>

namespace cues {
    struct SeqAln {
        std::string chr_name;
        int64_t ref_start;
        int64_t ref_end;
        int32_t query_start;
        int32_t query_end;
        bool is_reverse;
        uint8_t is_primary;
        uint8_t mapq;
        int64_t mlen;
        int64_t blen;
        uint32_t cigar_len;
        std::vector<uint32_t> cigar;

        inline io::SplitAln make_split() const {
            return io::SplitAln(std::string{}, chr_name, ref_start, ref_end, query_start, query_end, is_reverse, mapq);
        }

        inline std::string to_string() const {
            return std::format("{}@{}-{};{}-{}:{},{},{}", chr_name, ref_start, ref_end, query_start, query_end,
                               is_reverse, mapq, cigar_to_string());
        }

        inline std::string cigar_to_string() const {
            std::string cigar_str;
            for (auto i = 0; i < cigar_len; i++) {
                cigar_str += std::to_string(bam_cigar_oplen(cigar[i]));
                cigar_str += bam_cigar_opchr(bam_cigar_op(cigar[i]));
            }
            return cigar_str;
        }
    };

    class ReadAligner {
    public:
        explicit ReadAligner(const std::string &seq, const std::string &preset) {
            set_params(preset);
            const char* seqs[] = {seq.c_str()};
            index_ = mm_idx_str(index_opt_.w, index_opt_.k, index_opt_.flag&1, index_opt_.bucket_bits,
                                1, seqs, NULL);
            if (!index_) throw std::runtime_error("[MAP] Failed to build sequence index");
            mm_mapopt_update(&map_opt_, index_);
        }

        ~ReadAligner() {
            if (index_) mm_idx_destroy(index_);
        }
        ReadAligner(const ReadAligner&) = delete;
        ReadAligner& operator=(const ReadAligner&) = delete;

        std::vector<SeqAln> map(const std::string& query) {
            mm_tbuf_t* tbuf_ = mm_tbuf_init();
            std::vector<SeqAln> hits;
            int n_regs = 0;
            mm_reg1_t *regs = mm_map(index_, query.size(), query.c_str(), &n_regs, tbuf_, &map_opt_, "");
            for(auto i = 0; i < n_regs; i++) {
                const mm_reg1_t& hit = regs[i];
                if (!hit.p) continue;
                hits.push_back(SeqAln{
                        .chr_name =((index_->seq[hit.rid].name) ? std::string(index_->seq[hit.rid].name) : "local"),
                        .ref_start=hit.rs,
                        .ref_end=hit.re,
                        .query_start=hit.qs,
                        .query_end=hit.qe,
                        .is_reverse=static_cast<bool>(hit.rev),
                        .is_primary=(hit.id == hit.parent),
                        .mapq=static_cast<uint8_t>(hit.mapq),
                        .mlen=hit.mlen,
                        .blen=hit.blen,
                        .cigar_len=hit.p->n_cigar,
                        .cigar=std::vector<uint32_t>(hit.p->cigar, hit.p->cigar + hit.p->n_cigar)});
                free(hit.p);
            }
            free(regs);
            if (tbuf_) mm_tbuf_destroy(tbuf_);
            return hits;
        }

        private:
            ReadAligner() = default;
            void set_params(const std::string &preset) {
                mm_set_opt(NULL, &index_opt_, &map_opt_); // defaults
                if (!preset.empty() && preset != "custom" && preset != "global") { // standard preset
                    if (mm_set_opt(preset.c_str(), &index_opt_, &map_opt_) < 0) {
                        throw std::runtime_error("[MAP] Preset not supported: " + preset);
                    }
                }
                if (preset == "custom") { // TODO: adjust for "global"
                    mm_set_opt("sr", &index_opt_, &map_opt_);
                    index_opt_.w = 5;
                    index_opt_.k = 15;
                    map_opt_.min_cnt = 1;
                    map_opt_.min_dp_max = 10;
                    map_opt_.min_chain_score = 1;
                }
                map_opt_.flag |= MM_F_CIGAR;
                map_opt_.flag |= MM_F_EQX;
                map_opt_.best_n = 1;
                map_opt_.mid_occ = 1000;
            }
            mm_idxopt_t index_opt_;
            mm_mapopt_t map_opt_;
            mm_idx_t* index_;
    };

}
#endif //CUES_INCLUDE_CUES_ALIGN_H_
