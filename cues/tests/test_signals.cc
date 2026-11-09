/**
 * @brief SV signal extraction tests
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <array>
#include <filesystem>
#include <fstream>
#include <string>

#include "cues/io.h"
#include "cues/signals.h"

using namespace cues;
TEST_CASE("bam") {
    // BAM header
    const char* sam_header_str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    io::bam_header_t sam_hdr(sam_hdr_parse(strlen(sam_header_str), sam_header_str));
    // BAM file
    std::string bam_file_name = std::filesystem::temp_directory_path() / "tmpbamXXXXXX";
    io::BamWriter bam(bam_file_name, sam_hdr);
    std::vector<std::string> bam_entries = {
        "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r3\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r4\t0\tchr1\t101\t60\t5=10D5=\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r5\t0\tchr1\t101\t60\t5=10D5=\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r6\t0\tchr1\t101\t60\t5=5S\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r7\t0\tchr1\t101\t0\t5=5S\t*\t0\t0\tAAAAAAAAAA\t..........",
        "r8\t0\tchr1\t201\t60\t5=5I1=\t*\t0\t0\tAAAAAAAAAAA\t...........",
        "r9\t0\tchr1\t201\t60\t1S4=1I5S\t*\t0\t0\tAAAAAAAAAAA\t...........",
    };
    for (const auto& sam_record : bam_entries) {
        io::aln_ptr_t aln(bam_init1());
        kstring_t sam_entry = {0};
        kputs(sam_record.c_str(), &sam_entry);
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get());
        CHECK(ret >= 0);
        bam.write(aln);
    }
    bam.close();
    int ret = sam_index_build(bam_file_name.c_str(), 0);
    CHECK(ret >= 0);

    // Reference FASTA
    std::string fasta_file_name = std::filesystem::temp_directory_path() / "tmpfaXXXXXX";
    int fd = mkstemp(fasta_file_name.data());
    if (fd == -1) throw std::runtime_error("Error in temporary FASTA generation");
    close(fd);
    std::ofstream fasta_out(fasta_file_name);
    std::vector<std::pair<std::string, std::string>> ref_seq = {{"chr1",
                                                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG"
                                                                 + std::string(50, 'A')
                                                                 + std::string(400, 'C')}};
    for (const auto& [name, seq] : ref_seq) {
        fasta_out << ">" << name << "\n" << seq << "\n";
    }
    fasta_out.close();

    // Config
    cues::Cues::CuesConfig config;
    config.ref = fasta_file_name;
    config.min_mapq_frag = 0;
    config.run_recovery_subfrags = true;
    config.run_recovery_clipping = true;
    config.record_insertions = true;
    config.min_clipped_len = 30;
    config.min_discord_len = 30;
    config.record_inter_chromosomal_junctions = true;

    SUBCASE("cigar") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t5=100D5=\t*\t0\t0\tAAAAAAAAAA\t..........";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        CHECK(aln.get_ref_start() == 12344);
        CHECK(aln.get_ref_end() == 12453);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=12348, .pos2=12449,
                           .orientation=JOrientation::kLR,
                           .evidence_type=JSource::kCigar});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sr_ord12") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t10=10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t....................\t"
                                 "SA:Z:chr2,23456,+,10S10=,40,0;";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr2", .pos1=12353, .pos2=23455,
                                  .orientation=JOrientation::kLR,
                                  .evidence_type=JSource::kSplitRead});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sr_ord21") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t10=10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t....................\t"
                                 "SA:Z:chr1,12245,+,10S10=,40,0;";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=12244, .pos2=12353,
                                  .orientation=JOrientation::kRL,
                                  .evidence_type=JSource::kSplitRead});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sr_ord12_flip") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t10=10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t....................\t"
                                 "SA:Z:chr1,12445,-,10S10=,40,0;";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=12353, .pos2=12453,
                                  .orientation=JOrientation::kLL,
                                  .evidence_type=JSource::kSplitRead});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sf_ord21") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=50I5=\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=0, .pos2=50,
                                  .orientation=JOrientation::kRL,
                                  .evidence_type=JSource::kSubFragPair});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sf_flip_ins") {
        REQUIRE(sam_hdr != nullptr);
        std::string rec_seq = std::string("r1\t0\tchr1\t1\t60\t50=50X10=\t*\t0\t0\t") +
                      "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG" + std::string(50, 'T') + std::string(10, 'C') + "\t"
                      + "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG" + std::string(50, 'T') + std::string(10, 'C');
        const char *sam_record = rec_seq.c_str();
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=49, .pos2=99,
                                  .orientation=JOrientation::kLL,
                                  .evidence_type=JSource::kSubFragPair});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("sf_ord21_size_mismatch_small") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=52I5=\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAAAA\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
    }

    SUBCASE("sf_ord21_size_mismatch_large") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=60I5=\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGCCCCCCCCCCAAAAA\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGCCCCCCCCCCAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 0);
    }

    SUBCASE("sf_ord21_seq_mismatch_small") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=50I5=\t*\t0\t0\t"
                         "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTAATAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA\t"
                         "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
    }

    SUBCASE("sf_ord21_seq_mismatch_large") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=50I5=\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTAATAAAAAAATAACCCCCCGGCCGGGGGGGGCGAAAAA\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 0);
    }

    SUBCASE("ins_not_ord21") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=50I5=\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAATTTTGGAAAAAAAAAAAAAAAAAAAAAACCCGGGAAAAGGGGGGGGGGAAAAA\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAAAAA";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        CHECK(extractor.get_index().junctions.size() == 0);
        const auto& insertions = extractor.get_index().insertions;
        CHECK(insertions.size() == 1);
        CHECK(insertions[0].ref_pos == 50);
        CHECK(insertions[0].read_pos == 50);
    }

    SUBCASE("clip_ord21") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t1\t60\t50=50S\t*\t0\t0\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG\t"
                                 "AACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGAACCTTGGAATTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGG";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        io::BamAln aln("chr1");
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);

        cues::Cues cues(config);
        signals::SignalExtractor extractor(cues.get_config());
        extractor.process_aln(aln);
        const auto& junctions = extractor.get_index().junctions;
        CHECK(junctions.size() == 1);
        auto [junction, count] = *junctions.begin();
        Junction expected({.chr1 = "chr1", .chr2 = "chr1", .pos1=0, .pos2=49,
                                  .orientation=JOrientation::kRL,
                                  .evidence_type=JSource::kClip});
        CHECK(expected.to_string() == junction.to_string());
    }

    SUBCASE("stats_noop") {
        config.min_discord_len = 1;
        config.min_clipped_len = 1;
        cues::Cues cues(config);
        auto signal_index = cues.collect_bkps_from_bam(bam_file_name);
        std::vector<cues::RefSegment> segs = {cues::RefSegment("chr1", 0, 10), // fully covered, read length
                                                    cues::RefSegment("chr1", 5, 10), // fully covered
                                                    cues::RefSegment("chr1", 5, 15), // partial overlap
                                                    cues::RefSegment("chr1", 20, 30)}; // no overlap
        const auto stats = cues.collect_stats(bam_file_name, segs, {cues::AlnStatistic::kRD});
        CHECK(stats.size() == 4);
        CHECK(stats[0].at(cues::AlnStatistic::kRD) == 3);
        CHECK(stats[1].at(cues::AlnStatistic::kRD) == 3);
        CHECK(stats[2].at(cues::AlnStatistic::kRD) == 1.5);
        CHECK(stats[3].at(cues::AlnStatistic::kRD) == 0);
    }

    SUBCASE("stats_sv") {
        config.min_discord_len = 1;
        config.min_clipped_len = 1;
        cues::Cues cues(config);
        auto signal_index = cues.collect_bkps_from_bam(bam_file_name);
        std::vector<cues::RefSegment> segs = {cues::RefSegment("chr1", 100, 105), // fully covered
                                                    cues::RefSegment("chr1", 100, 110), // internal deletion
                                                    cues::RefSegment("chr1", 100, 125), // partial overlap
                                                    cues::RefSegment("chr1", 115, 120),
                                                    cues::RefSegment("chr1", 200, 210), // insertions
                                                    cues::RefSegment("chr1", 195, 205)}; // insertions boundary
        const auto stats = cues.collect_stats(bam_file_name, segs,
        {cues::AlnStatistic::kRD, cues::AlnStatistic::kRDAll, cues::AlnStatistic::kClipRight,
        cues::AlnStatistic::kClipLeft, cues::AlnStatistic::kIns});
        CHECK(stats[0].at(cues::AlnStatistic::kRD) == 3);
        CHECK(stats[0].at(cues::AlnStatistic::kRDAll) == 4);
        CHECK(stats[1].at(cues::AlnStatistic::kRD) == 1.5);
        CHECK(stats[1].at(cues::AlnStatistic::kClipRight) == 1);
        CHECK(stats[2].at(cues::AlnStatistic::kRD) == 1);
        CHECK(stats[3].at(cues::AlnStatistic::kRD) == 2);
        CHECK(stats[4].at(cues::AlnStatistic::kIns) == 2);
        CHECK(stats[4].at(cues::AlnStatistic::kClipLeft) == 1);
        CHECK(stats[4].at(cues::AlnStatistic::kClipRight) == 1);
        CHECK(stats[5].at(cues::AlnStatistic::kIns) == 1);
    }
}