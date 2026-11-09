// Copyright (c) 2026 Victoria Popic
// AGPL-3.0-or-later

#include "cues/cues.h"
#include "cues/signals.h"
#include "cues/io.h"
#include <sys/resource.h>
#include <algorithm>
#include <cmath>
#include <omp.h>

namespace cues {
    static void check_bam(const io::aln_buffer_t& batch) {
        for (const auto& aln : batch) {
            if ((aln.get_flag() & (BAM_FUNMAP | BAM_FSECONDARY)) == 0 && aln.is_hard_clipped()) {
                throw std::runtime_error("[BAM error] cues requires soft-clipped BAM");
            }
        }
    }
    BreakpointStore Cues::collect_bkps_from_bam(const std::string &bam_file_name) {
        CUES_LOG("Collecting candidate junctions from " << bam_file_name);
        io::BamReader input_bam(bam_file_name, config_.aln_batch_size, config_.max_alns_to_process);
        std::vector<signals::SignalExtractor> index_per_thr;
        index_per_thr.reserve(config_.n_threads);
        std::unordered_map<std::string, ReadAligner> chr2aligner;

        if (config_.run_recovery_ins || config_.run_recovery_missing_segments) {
            io::Fasta fasta(config_.ref);
            for (const auto& chr_name : fasta.get_chr_names()) {
                if (signals::SignalExtractor::skip_chr(chr_name)) continue;
                chr2aligner.try_emplace(chr_name, fasta.fetch(chr_name), "sr");
            }
        }
        for (auto i = 0; i < config_.n_threads; i++) index_per_thr.emplace_back(config_, chr2aligner);
        while (true) {
            auto &batch = input_bam.read_next_batch();
            if (batch.empty()) break;
            check_bam(batch);
            const auto aln_per_thread = static_cast<int>(std::ceil(static_cast<double>(batch.size()) / config_.n_threads));
            #pragma omp parallel for default(none) shared(batch, aln_per_thread, index_per_thr, input_bam, std::cout)
            for (size_t worker_start = 0; worker_start < batch.size(); worker_start += aln_per_thread) {
                const auto worker_end = std::min(worker_start + aln_per_thread, batch.size());
                const auto tid = omp_get_thread_num();
                for (auto aln_id = worker_start; aln_id < worker_end; aln_id++) {
                    index_per_thr[tid].process_aln(batch[aln_id]);
                }
            }
            CUES_LOG("Processed " << input_bam.get_n_processed() << " reads");
        }
        // merge per-thread outputs
        signals::SignalIndex combined_index;
        BreakpointStore output;
        for (auto& index : index_per_thr) {
            for (const auto& [j, count] : index.get_index().junctions) {
                combined_index.junctions[j] += count;
            }
            output.insertions.insert(output.insertions.end(), index.get_index().insertions.begin(),
                                     index.get_index().insertions.end());
            output.aln_patches.insert(output.aln_patches.end(), index.get_index().aln_patches.begin(),
                                     index.get_index().aln_patches.end());
        }

        for (const auto& [j, count] : combined_index.junctions) {
            Junction j_out = j;
            j_out.support = count;
            output.junctions.push_back(j_out);
        }
        CUES_LOG("Extracted " << combined_index.junctions.size() << " junctions");
        CUES_LOG("Total number of reads processed: " << input_bam.get_n_processed());
        return output;
    }

    std::vector<aln_stats_map_t> Cues::collect_stats(const std::string& bam_file_name,
                                                     const std::vector<RefSegment>& segments,
                                                     const std::unordered_set<AlnStatistic>& stats_to_compute,
                                                     const std::vector<AlnPatch>& aln_patches) {
        omp_set_num_threads(config_.n_threads);
        std::vector<aln_stats_map_t> stats(segments.size());
        signals::AlnStatsGenerator stats_generator(config_, aln_patches);
        const auto seg_per_thr = static_cast<int>(std::ceil(static_cast<double>(segments.size()) / config_.n_threads));
        #pragma omp parallel default(none) \
            shared(bam_file_name, segments, stats, stats_generator, seg_per_thr, stats_to_compute, std::cout)
        {
            io::BamReader input_bam(bam_file_name);
            #pragma omp for schedule(static, 1)
            for (size_t seg_id = 0; seg_id < segments.size(); ++seg_id) {
                const auto& segment = segments[seg_id];
                io::aln_buffer_t alns = input_bam.read_region_overlap(segment.chr, segment.start, segment.end,
                                                                      config_.max_seg_alns_to_fetch);
                stats[seg_id] = stats_generator.compute_stats(alns, segment, stats_to_compute);
            }
        }
        return stats;
    }
}
