#ifndef CUES_INCLUDE_CUES_SIGNALS_H_
#define CUES_INCLUDE_CUES_SIGNALS_H_

#include "cues/align.h"
#include "cues/cues.h"
#include "cues/io.h"
#include <algorithm>
#include <format>
#include <stdlib.h>
#include <unordered_map>
#include <unordered_set>

namespace cues {
    namespace signals {
        struct SignalIndex {
            junctions_map_t junctions;
            int n_junctions;
            std::vector<NovelSeq> insertions;
            std::unordered_map<std::string, std::vector<Junction>> chr2cdels; // sorted by reference position
            std::unordered_map<std::string, std::vector<NovelSeq>> chr2cins;  // sorted by reference position
            std::vector<AlnPatch> aln_patches;

            inline void extract_cdels(BreakpointStore& sv_cues) {
                for (const auto& j : sv_cues.junctions) {
                    if (j.evidence_type == JSource::kCigar) chr2cdels[j.chr1].push_back(j);
                }
                for (auto& [chr, cdels] : chr2cdels) {
                    std::sort(cdels.begin(), cdels.end(), [](const auto& j1, const auto& j2) { return j1.pos1 < j2.pos1; });
                }
            }

            inline void extract_cins(BreakpointStore& sv_cues) {
                for (const auto& i : sv_cues.insertions) {
                    chr2cins[i.chr].push_back(i);
                }
                for (auto& [chr, cins] : chr2cins) {
                    std::sort(cins.begin(), cins.end(), [](const auto& i1, const auto& i2) { return i1.ref_pos < i2.ref_pos; });
                }
            }
        };

        class SignalExtractor {
        public:
            explicit SignalExtractor(Cues::CuesConfig config):
            config_(config), fasta_(config.ref), chr2aligner_(nullptr) {}
            explicit SignalExtractor(Cues::CuesConfig config, std::unordered_map<std::string, ReadAligner>& chr2aligner):
            config_(config), fasta_(config.ref), chr2aligner_(&chr2aligner) {}
            void process_aln(io::BamAln& aln);
            void process_read(io::FastqRead& read);
            inline SignalIndex &get_index() { return index_; }
            inline std::unordered_set<std::string> &get_reads() { return seen_reads_; }
            static inline bool skip_chr(const std::string& chr_name) {
                const std::string chr_type = chr_name.starts_with("chr") ? chr_name.substr(3) : chr_name;
                if (chr_type == "X" || chr_type == "Y") return false;
                return chr_type.empty() || !std::all_of(chr_type.begin(), chr_type.end(), ::isdigit);
            }
        private:
	        bool skip_aln(const io::BamAln& aln) const;
            inline bool record_junction(const Junction& j) {
                if (skip_chr(j.chr1) || skip_chr(j.chr2)) {
                    return false;
                }
                if (j.chr1 == j.chr2) {
		            if (j.len() >= config_.min_discord_len && j.len() <= config_.max_discord_len) {
                        index_.junctions[j] += 1;
                        index_.n_junctions++;
                        return true;
                    }
                } else if (config_.record_inter_chromosomal_junctions) {
                    index_.junctions[j] += 1;
                    index_.n_junctions++;
                    return true;
                }
                return false;
            }
            inline void record_insertion(const NovelSeq& i) {
                if (!config_.record_insertions) return;
                if (i.len >= config_.min_discord_len) index_.insertions.push_back(i);
            }
            inline void record_aln_patch(const AlnPatch& p) {
                if (!config_.record_aln_patches) return;
                index_.aln_patches.push_back(p);
            }
            bool process_cigar(io::BamAln& aln);
            void process_split_read(io::BamAln& aln);
            void process_clipping(io::BamAln& aln);
            void process_ins(io::BamAln& aln);
            std::vector<SeqAln> realign(const std::string& chr_name, const std::string& target_seq);
            void run_signal_recovery(io::BamAln& aln);
            Junction make_frag_junction(const SeqAln &aln1, const SeqAln &aln2, const io::BamAln& parent_aln) const;
            static Junction make_sr_junction(const io::SplitAln& aln1, const io::SplitAln& aln2);
            bool verify_and_refine_junction(Junction& junction, io::BamAln& parent_aln);
			void trim_query_overlaps(std::vector<io::SplitAln>& split_alns);
            Cues::CuesConfig config_;
            io::Fasta fasta_;
            SignalIndex index_;
            std::unordered_set<std::string> seen_reads_;
            std::unordered_map<std::string, ReadAligner>* chr2aligner_;
        };

        class AlnStatsGenerator {
        public:
            explicit AlnStatsGenerator(Cues::CuesConfig config) : config_(config) {}
 			AlnStatsGenerator(Cues::CuesConfig config, const std::vector<AlnPatch>& patches)
			: config_(config) {
                for (const auto& p : patches) patch_index_.insert_or_assign(p.aln_id, p);
			}
            aln_stats_map_t compute_stats(io::aln_buffer_t &alns, const RefSegment &segment,
                                          const std::unordered_set <AlnStatistic> &stats_to_compute) const;

        private:
            Cues::CuesConfig config_;
            std::unordered_map<std::string, AlnPatch> patch_index_;
        };
    }
}
#endif
