// Copyright (c) 2026 Victoria Popic
// AGPL-3.0-or-later

#include "cues/io.h"
#include "cues/signals.h"
#include <algorithm>
#include <numeric>

namespace cues {
    namespace signals {
        bool SignalExtractor::skip_aln(const io::BamAln& aln) const {
            std::string chr_name = aln.get_chr_name();
            if (skip_chr(chr_name)) return true;
            const auto flag = aln.get_flag();
            if ((flag & BAM_FUNMAP) || (flag & BAM_FSECONDARY) || (aln.get_mapq() < config_.min_mapq)) return true;
            return aln.get_discordance() > config_.max_discordance;
        }

        void SignalExtractor::process_aln(io::BamAln& aln) {
            if (skip_aln(aln)) return;
            const bool is_divergent = process_cigar(aln);
            if (aln.is_split()) {
                process_split_read(aln);
            } else if (config_.run_recovery_clipping and aln.is_clipped()) {
                process_clipping(aln);
            }
            if (!is_divergent and aln.cins.empty()) return;
            if (config_.run_recovery_subfrags) run_signal_recovery(aln);
            if (config_.run_recovery_ins and !aln.cins.empty()) process_ins(aln);
        }

        bool SignalExtractor::process_cigar(io::BamAln& aln) {
            const uint32_t *cigar = aln.get_cigar();
            ref_pos_t ref_pos = aln.get_ref_start();
            int read_pos = 0;
            int discordant_window_len = 0;
            bool is_divergent = false;
            for (hts_pos_t i = 0; i < aln.get_cigar_len(); ++i) {
                const int op = bam_cigar_op(cigar[i]);
                const int op_len = bam_cigar_oplen(cigar[i]);
                switch (op) {
                    case BAM_CINS:
                        if (op_len >= config_.min_discord_len) {
                            aln.cins.emplace_back(ref_pos, read_pos, op_len, false);
                        }
                        record_insertion(NovelSeq({.chr=aln.get_chr_name(), .ref_pos=ref_pos,
                                                   .read_pos=read_pos, .len=op_len}));
                        if (op_len < config_.min_discord_len) {
                            discordant_window_len += op_len;
                        }
                        read_pos += op_len;
                        break;
                    case BAM_CDEL:
                        record_junction(Junction(
                                {.chr1 = aln.get_chr_name(),
                                 .chr2 = aln.get_chr_name(),
                                 .pos1=ref_pos-1,
                                 .pos2=ref_pos + op_len,
                                 .orientation=JOrientation::kLR,
                                 .evidence_type=JSource::kCigar,
                                 }));
                        if (op_len < config_.min_discord_len) {
                            discordant_window_len += op_len;
                        }
                        ref_pos += op_len;
                        break;
                    case BAM_CDIFF:
                        discordant_window_len += op_len;
                        ref_pos += op_len;
                        read_pos += op_len;
                        break;
                    case BAM_CMATCH:
                    case BAM_CEQUAL:
                        if (op_len > config_.min_accord_len) {
                            discordant_window_len = 0;
                        }
                        ref_pos += op_len;
                        read_pos += op_len;
                        break;
                    case BAM_CHARD_CLIP:
                        break;
                    case BAM_CSOFT_CLIP:
                        read_pos += op_len;
                        break;
                    default:
                        break;
                }
                if (discordant_window_len >= config_.min_discord_len) {
                    is_divergent = true;
                }
            }
            return is_divergent;
        }

        void SignalExtractor::process_split_read(io::BamAln& aln) {
            std::vector<io::SplitAln> split_alns;
            split_alns.push_back(aln.make_split());
            for (auto& sa : aln.get_supplementary_alns()) {
                if (sa.mapq == 0) continue; // TODO
                split_alns.push_back(sa);
            }
            std::ranges::sort(split_alns, {}, &io::SplitAln::query_start);
            if (split_alns[0].query_start != aln.get_query_start()) return;

            if (config_.run_recovery_missing_segments) {
                const auto read_length = aln.get_seq().size();
                std::vector<io::SplitAln> recovered_segments;
                for (int i = 1; i < split_alns.size(); i++) {
                    const auto segment_length = split_alns[i].query_start - 1 - split_alns[i - 1].query_end;
                    if (segment_length >= config_.min_discord_len) {
                        const auto segment_start = !aln.is_reverse() ? split_alns[i - 1].query_end + 1 : read_length - split_alns[i].query_start;
                        auto segment_seq = aln.get_seq().substr(segment_start, segment_length);
                        if (aln.is_reverse()) {
                            cues::io::revcomp(segment_seq, aln.get_name());
                        }
                        const auto hits = realign(aln.get_chr_name(), segment_seq);
                        if (hits.size() != 1) continue;
                        auto hit = hits[0].make_split();
                        hit.chr_name = aln.get_chr_name();
                        hit.query_start = split_alns[i - 1].query_end + 1;
                        hit.query_end = hit.query_start + segment_length - 1;
                        recovered_segments.push_back(hit);
                    }
                }
                split_alns.insert(split_alns.end(), recovered_segments.begin(), recovered_segments.end());
                std::ranges::sort(split_alns, {}, &io::SplitAln::query_start);
            }

            trim_query_overlaps(split_alns);

            for (int i = 1; i < split_alns.size(); i++) {
                // do not create a junction when a read segment is missing
                if (split_alns[i].query_start > split_alns[i-1].query_end + config_.max_segment_missing) {
                    continue;
                }
                const auto& junction = make_sr_junction(split_alns[i-1], split_alns[i]);
                record_junction(junction);
            }
        }

        void SignalExtractor::run_signal_recovery(io::BamAln& aln) {
            const std::string ref_seq = fasta_.fetch(aln.get_chr_name(), aln.get_ref_start(), aln.get_ref_end());
            auto aligner = ReadAligner(ref_seq, "custom");
            std::optional<SeqAln> last_hit;
            for (size_t i = 0; i + config_.frag_len + 1 <= aln.get_seq().size(); i += config_.frag_len) {
                std::string fragment = aln.get_seq().substr(i, config_.frag_len);
                auto hits = aligner.map(fragment);
                if (hits.empty()) {
                    last_hit = std::nullopt;
                    continue;
                }
                const SeqAln &hit = hits.front();
                if (hit.mapq < config_.min_mapq_frag) {
                    last_hit = std::nullopt;
                    continue;
                }
                if (!last_hit.has_value()) {
                    last_hit = hit;
                    continue;
                }
                auto junction = make_frag_junction(last_hit.value(), hit, aln);
                if (junction.len() >= config_.min_discord_len && junction.len() <= config_.max_discord_len) {
                    if (verify_and_refine_junction(junction, aln)) {
                        record_junction(junction);
                    }
                }
                last_hit = hit;
            }
        }

        Junction SignalExtractor::make_frag_junction(const SeqAln &aln1, const SeqAln &aln2,
                                                 const io::BamAln& parent_aln) const {
            const ref_pos_t ref_start1 = parent_aln.get_ref_start() + aln1.ref_start;
            const ref_pos_t ref_end1 = parent_aln.get_ref_start() + aln1.ref_end;
            const ref_pos_t ref_start2 = parent_aln.get_ref_start() + aln2.ref_start;
            const ref_pos_t ref_end2 = parent_aln.get_ref_start() + aln2.ref_end;
            ref_pos_t j1_shift = 0;
            ref_pos_t j2_shift = 0;
            ref_pos_t junction_start;
            if (!aln1.is_reverse) {
                junction_start = ref_end1-1;
                if (aln1.query_end < config_.frag_len) j2_shift = config_.frag_len - aln1.query_end;
                if (aln2.query_start > 0) j1_shift = aln2.query_start;
            } else {
                junction_start = ref_start1;
                if (aln1.query_end < config_.frag_len) j2_shift = -(config_.frag_len - aln1.query_end);
                if (aln2.query_start > 0) j1_shift = -aln2.query_start;
            }
            junction_start += j1_shift;
            const ref_pos_t junction_end = aln2.is_reverse ? ref_end2 - 1 : ref_start2 + j2_shift;
            const bool reorder = junction_start > junction_end;
            return Junction({
                    .chr1 = parent_aln.get_chr_name(),
                    .chr2 = parent_aln.get_chr_name(),
                    .pos1 = std::min(junction_start, junction_end),
                    .pos2 = std::max(junction_start, junction_end),
                    .orientation = aln1.is_reverse == aln2.is_reverse ?
                                   (((!reorder and !aln1.is_reverse) || (reorder and aln1.is_reverse)) ?
                                   JOrientation::kLR : JOrientation::kRL)
                                   : ((junction_start > aln1.ref_start + parent_aln.get_ref_start()) ? JOrientation::kLL : JOrientation::kRR),
                    .evidence_type = JSource::kSubFragPair
            });
        }

        Junction SignalExtractor::make_sr_junction(const io::SplitAln& aln1, const io::SplitAln& aln2) {
            const ref_pos_t junction_start = aln1.is_reverse ? aln1.ref_start : aln1.ref_end;
            const ref_pos_t junction_end = aln2.is_reverse ? aln2.ref_end : aln2.ref_start;
            const bool reorder = junction_start > junction_end;
            return Junction({.chr1 = reorder ? aln2.chr_name : aln1.chr_name,
                             .chr2 = reorder ? aln1.chr_name : aln2.chr_name,
                             .pos1 = reorder ? junction_end : junction_start,
                             .pos2 = reorder ? junction_start : junction_end,
                             .orientation = aln1.is_reverse == aln2.is_reverse ?
                                     (((!reorder and !aln1.is_reverse) || (reorder and aln1.is_reverse)) ?
                                     JOrientation::kLR : JOrientation::kRL)
                                     : ((junction_start > aln1.ref_start) ? JOrientation::kLL : JOrientation::kRR),
                             .evidence_type = JSource::kSplitRead
                            });
        }

        void SignalExtractor::trim_query_overlaps(std::vector<io::SplitAln>& split_alns) {
            std::vector<std::optional<AlnPatch>> patches(split_alns.size());
            for (auto i = 1; i < split_alns.size(); i++) {
                auto& first  = split_alns[i - 1];
                auto& second = split_alns[i];
                const auto overlap = static_cast<int64_t>(first.query_end)
                                   - static_cast<int64_t>(second.query_start);
                if (overlap <= 0) continue;
                const auto junc1_ref = first.is_reverse  ? first.ref_start : first.ref_end;
                const auto junc2_ref = second.is_reverse ? second.ref_end  : second.ref_start;
                const bool trim_first = first.is_reverse ? (junc1_ref < junc2_ref) : (junc1_ref > junc2_ref);
                const auto target_idx = trim_first ? i - 1 : i;
                auto& target = split_alns[target_idx];
                const auto len = static_cast<int64_t>(target.query_end)
                               - static_cast<int64_t>(target.query_start);
                if (overlap >= len) continue;
                const bool trim_ref_end = trim_first ? !target.is_reverse : target.is_reverse;
                if (!patches[target_idx].has_value()) {
                    patches[target_idx] = AlnPatch{io::make_aln_uid(target), RefSegment{target.chr_name, 0, 0}};
                }
                if (trim_ref_end) target.ref_end -= overlap;
                else target.ref_start += overlap;
                if (trim_first) target.query_end -= overlap;
                else target.query_start += overlap;
                patches[target_idx]->new_span.start = target.ref_start;
                patches[target_idx]->new_span.end   = target.ref_end;
            }
            for (const auto& p : patches) {
                if (p.has_value()) record_aln_patch(*p);
            }
        }

        bool SignalExtractor::verify_and_refine_junction(Junction& junction, io::BamAln& aln) {
            double size_sim = 0.7;
            double seq_sim = 0.9;
            int max_custom_len = 500;
            if ((junction.orientation == JOrientation::kLL) || (junction.orientation == JOrientation::kRR)) {
                std::string ref_seq = fasta_.fetch(aln.get_chr_name(), junction.pos1, junction.pos2);
                auto aligner = ReadAligner(aln.get_seq(), ref_seq.size() <= max_custom_len ? "custom" : "sr");
                auto hits = aligner.map(ref_seq);
                if (hits.empty() or !hits[0].is_reverse) return false;
                const SeqAln &hit = hits.front();
                const double size_ratio = double(std::min(hit.blen, junction.len()))/std::max(hit.blen, junction.len());
                const double seq_ratio = double(hit.mlen) / hit.blen;
                if (size_ratio < size_sim || seq_ratio < seq_sim) return false;
                return true;
            } else if (junction.orientation == JOrientation::kRL) {
                std::vector<std::tuple<double, ref_pos_t, ref_pos_t, ref_pos_t, size_t>> insertion_match;
                const double slack = (1.0 - size_sim) * double(junction.len());
                for (auto i = 0; i < aln.cins.size(); i++) {
                    const auto& [ref_pos, read_pos, len, matched] = aln.cins[i];
                    if (double(ref_pos) < double(junction.pos1) - slack or
                        double(ref_pos) > double(junction.pos2) + slack) continue;
                    const double size_ratio = double(std::min(len, junction.len()))/std::max(len, junction.len());
                    if (size_ratio < size_sim) continue;
                    insertion_match.emplace_back(size_ratio, ref_pos, read_pos, len, i);
                }
                if (insertion_match.empty()) return false;

                const auto& [_, ref_pos, read_pos, len, idx] = *std::max_element(
                        insertion_match.begin(), insertion_match.end(), [](const auto& a, const auto& b)
                        {return get<0>(a) < get<0>(b);});
                const ref_pos_t ref_offset = std::max(ref_pos_t(0), ref_pos - len);
                std::string ref_seq = fasta_.fetch(aln.get_chr_name(), ref_offset,
                                                   std::min(ref_pos + len, fasta_.get_chr_len(aln.get_chr_name())));
                auto aligner = ReadAligner(ref_seq, len <= max_custom_len ? "custom" : "sr");
                ref_pos_t mlen = 0;
                std::vector<ref_pos_t> coords;
                auto hits = aligner.map(aln.get_seq().substr(read_pos, len));
                for (const auto& hit : hits) {
                    if (!hit.is_primary) continue;
                    mlen += hit.mlen;
                    coords.push_back(hit.ref_start);
                    coords.push_back(hit.ref_end);
                }
                if (aln_stats_value_t(mlen)/aln_stats_value_t(len) < seq_sim) return false;
                junction.pos1 = ref_offset + std::ranges::min(coords);
                junction.pos2 = ref_offset + std::ranges::max(coords);
                std::get<3>(aln.cins[idx]) = true; // mark the matching insertion in the index
                return true;
            }
            return false;
        }

        void SignalExtractor::process_clipping(io::BamAln& aln) {
            double seq_sim = 0.9;
            auto aligner = ReadAligner(fasta_.fetch(aln.get_chr_name(), aln.get_ref_start(), aln.get_ref_end()), "sr");
            const auto& read_seq = aln.get_seq();
            for (std::string side : {"left", "right"}) {
                const auto clip_len = aln.get_clipping(side);
                if (clip_len < config_.min_discord_len) continue;
                auto hits = aligner.map(side == "left" ? read_seq.substr(0, clip_len) : read_seq.substr(read_seq.size() - clip_len));
                if (hits.empty() || double(hits[0].mlen)/clip_len < seq_sim) continue;
                const auto&hit = hits.front();
                if (hit.is_reverse != aln.is_reverse()) continue;
                ref_pos_t junction_start, junction_end;
                if (side == "left") {
                    junction_start = aln.get_ref_start() + (hit.is_reverse ? hit.ref_start : hit.ref_end);
                    junction_end = !aln.is_reverse() ? aln.get_ref_start() : aln.get_ref_end();
                } else {
                    junction_start = !aln.is_reverse() ? aln.get_ref_end() : aln.get_ref_start();
                    junction_end = aln.get_ref_start() + (hit.is_reverse ? hit.ref_end : hit.ref_start);

                }
                record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                        .pos1 = std::min(junction_start, junction_end),
                        .pos2 = std::max(junction_start, junction_end),
                        .orientation = JOrientation::kRL,
                        .evidence_type = JSource::kClip});
            }
        }

        void SignalExtractor::process_ins(io::BamAln& aln) {
            for (auto i = 0; i < aln.cins.size(); i++) {
                const auto &[ref_pos, read_pos, len, matched] = aln.cins[i];
                if (matched) continue;
                const std::string read_seq = aln.get_seq().substr(read_pos, len);
                auto hits = realign(aln.get_chr_name(), read_seq);
                if (hits.size() != 1) continue;
                const auto hit = hits[0];
                std::get<3>(aln.cins[i]) = true;
                if (!hit.is_reverse) {
                    if (hit.ref_start < ref_pos) {
                        record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                                .pos1 = hit.ref_start,
                                .pos2 = ref_pos - 1,
                                .orientation = JOrientation::kRL,
                                .evidence_type = JSource::kIns});
                        record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                                .pos1 = std::min(hit.ref_end, ref_pos),
                                .pos2 = std::max(hit.ref_end, ref_pos),
                                .orientation = JOrientation::kLR,
                                .evidence_type = JSource::kIns});
                    } else {
                        record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                                .pos1 = std::min(ref_pos, hit.ref_end),
                                .pos2 = std::max(ref_pos, hit.ref_end),
                                .orientation = JOrientation::kRL,
                                .evidence_type = JSource::kIns});
                        record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                                .pos1 = std::min(ref_pos - 1, hit.ref_start),
                                .pos2 = std::max(ref_pos - 1, hit.ref_start),
                                .orientation = JOrientation::kLR,
                                .evidence_type = JSource::kIns});
                    }
                } else {
                    record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                            .pos1 = std::min(hit.ref_end, ref_pos - 1),
                            .pos2 = std::max(hit.ref_end, ref_pos - 1),
                            .orientation = JOrientation::kLL,
                            .evidence_type = JSource::kIns});
                    record_junction(Junction{.chr1 = aln.get_chr_name(), .chr2 = aln.get_chr_name(),
                            .pos1 = std::min(hit.ref_start, ref_pos),
                            .pos2 = std::max(hit.ref_start, ref_pos),
                            .orientation = JOrientation::kRR,
                            .evidence_type = JSource::kIns});
                }
            }
        }

        std::vector<SeqAln> SignalExtractor::realign(const std::string& chr_name, const std::string& target_seq) {
            double seq_sim = 0.95;
            double size_sim = 0.95;
            auto &aligner = chr2aligner_->at(chr_name);
            std::vector<SeqAln> hits;
            int n_hits_processed = 0;
            for (const auto &hit: aligner.map(target_seq)) { // find all the possible matches
                n_hits_processed++;
                if (n_hits_processed > 10) break;
                if (!hit.is_primary) continue;
                if (hit.mapq < config_.min_mapq_frag) continue;
                const ref_pos_t len = target_seq.size();
                const double size_ratio = double(std::min(hit.blen, len)) / double(std::max(hit.blen, len));
                const double seq_ratio = double(hit.mlen) / double(hit.blen);
                if (size_ratio < size_sim || seq_ratio < seq_sim) continue;
                if (!hits.empty()) break;
                hits.push_back(hit);
            }
            return hits;
        }

        aln_stats_map_t AlnStatsGenerator::compute_stats(io::aln_buffer_t& alns, const RefSegment& segment,
                                                         const std::unordered_set<AlnStatistic>& stats_to_compute) const {
            // TODO: error on zero-len segments
            aln_stats_map_t stats;
            for (auto& stat : stats_to_compute) {
                stats[stat] = 0;
            }

			size_t n_alns = 0;
            for (auto& aln : alns) {
                // apply alignment trimming fixes
                ref_pos_t ref_start = aln.get_ref_start();
                ref_pos_t ref_end   = aln.get_ref_end();
                const auto it = patch_index_.find(io::make_aln_uid(aln));
                if (it != patch_index_.end()) {
                    ref_start = it->second.new_span.start;
                    ref_end   = it->second.new_span.end;
                }
				if (ref_end < segment.start || ref_start >= segment.end) continue;

                // find the length of the alignment and segment overlap
                const auto aln_overlap = std::min(ref_end + 1, segment.end) - std::max(ref_start, segment.start);
                // find the length of cigar deletions, insertions, mismatches overlapping the segment
                auto [cins_total, cdel_total, cdiff_total] = (it != patch_index_.end())
                ? aln.parse_cigar_in_range(std::max(ref_start, segment.start),
                                           std::min(ref_end + 1, segment.end))
                        : aln.parse_cigar_in_range(segment.start, segment.end);

                // read depth (all)
                if (stats_to_compute.contains(AlnStatistic::kRDAll)) {
                    stats[AlnStatistic::kRDAll] += (aln_stats_value_t) aln_overlap;
                    stats[AlnStatistic::kRDAll] -= (aln_stats_value_t) cdel_total;

                }
                if (aln.get_mapq() < config_.min_mapq) continue;
				n_alns++;

                // read depth (high confidence)
                if (stats_to_compute.contains(AlnStatistic::kRD)) {
                    stats[AlnStatistic::kRD] += (aln_stats_value_t) aln_overlap;
                    stats[AlnStatistic::kRD] -= (aln_stats_value_t) cdel_total;
                }

                // insertion presence
                if (stats_to_compute.contains(AlnStatistic::kIns)) {
                    if (cins_total >= config_.min_discord_len) {
                        stats[AlnStatistic::kIns] += 1;
                    }
                }

                if (stats_to_compute.contains(AlnStatistic::kDiv)) {
                    stats[AlnStatistic::kDiv] += (cins_total + cdel_total + cdiff_total) / (aln_stats_value_t) aln_overlap;
                }

                // clipping
                // count it only if it occurs in this segment and has the required length
                if (stats_to_compute.contains(AlnStatistic::kClipLeft)) {
                    if (aln.get_clipping("left") >= config_.min_clipped_len) {
                        if (segment.start <= ref_start && ref_start < segment.end) {
                            stats[AlnStatistic::kClipLeft] += 1;
                        }
                    }
                }
                if (stats_to_compute.contains(AlnStatistic::kClipRight)) {
                    if (aln.get_clipping("right") >= config_.min_clipped_len) {
                        if (segment.start <= ref_end && ref_end < segment.end) {
                            stats[AlnStatistic::kClipRight] += 1;
                        }
                    }
                }
            }

            // update depth to average over the segment
            if (stats_to_compute.contains(AlnStatistic::kRDAll)) {
                stats[AlnStatistic::kRDAll] /= (aln_stats_value_t) segment.len();
            }

            if (stats_to_compute.contains(AlnStatistic::kRD)) {
                stats[AlnStatistic::kRD] /= (aln_stats_value_t) segment.len();
            }

            if (stats_to_compute.contains(AlnStatistic::kDiv) and n_alns > 0) {
                stats[AlnStatistic::kDiv] /= (aln_stats_value_t) n_alns;
            }
            return stats;
        }
    }
}
