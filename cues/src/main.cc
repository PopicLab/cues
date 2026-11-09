// Copyright (c) 2026 Victoria Popic
// AGPL-3.0-or-later

#include <filesystem>
#include <iostream>
#include <omp.h>
#include <string>

#include "cues/cues.h"

int main(int argc, char const ** argv) {
    std::cout << "********\n";
    std::cout << "* cues *\n";
    std::cout << "********\n";
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <ref_file_path> <bam_file_path>\n";
        return EXIT_FAILURE;
    }
    std::string ref_file_path = argv[1];
    std::string bam_file_path = argv[2];
    if (!std::filesystem::exists(ref_file_path)) {
        std::cerr << "FASTA file " << ref_file_path << " does not exist\n";
        return EXIT_FAILURE;
    }
    if (!std::filesystem::exists(bam_file_path)) {
        std::cerr << "BAM file " << bam_file_path << " does not exist\n";
        return EXIT_FAILURE;
    }

    const double start_time = omp_get_wtime();
    try {
        cues::Cues::CuesConfig config;
        config.ref = ref_file_path;
        config.max_alns_to_process = 1000000;
        config.aln_batch_size = 1000000;
        config.n_threads = 24;
        config.min_discord_len = 40;
        config.record_aln_patches = true;
        cues::Cues cues(config);
        auto index = cues.collect_bkps_from_bam(bam_file_path);
        std::cout << "Total adj time: " << omp_get_wtime() - start_time << " seconds\n";
        std::vector<cues::RefSegment> segs;
        const int radius = 20;
        for (auto &j : index.junctions) {
            segs.emplace_back(j.chr1, j.pos1-radius, j.pos1);
            segs.emplace_back(j.chr2, j.pos2, j.pos2+radius);
        }
        const auto signals = cues.collect_stats(bam_file_path, segs, {cues::AlnStatistic::kRD,
            cues::AlnStatistic::kClipRight, cues::AlnStatistic::kClipLeft, cues::AlnStatistic::kIns}, index.aln_patches);
        std::cout << "Total stats time: " << omp_get_wtime() - start_time << " seconds\n";
        for (int i = 0; i < signals.size(); i++) {
            std::cout << segs[i].chr << "[" << segs[i].start << "-" << segs[i].end
            << "] RD=" << signals[i].at(cues::AlnStatistic::kRD) <<"\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "[ERROR] " << e.what() << std::endl;
        std::exit(EXIT_FAILURE);
    }
    return EXIT_SUCCESS;
}
