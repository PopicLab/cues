/**
 * @brief BAM IO tests (alignment representation, retrieval, parsing, inference)
 */

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include <doctest/doctest.h>
#include <array>
#include <string>

#include "cues/io.h"

using namespace cues::io;
TEST_CASE("bam_io") {
    const char* sam_header_str = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:248956422\n";
    bam_header_t sam_hdr(sam_hdr_parse(strlen(sam_header_str), sam_header_str));

    SUBCASE("basic") {
        REQUIRE(sam_hdr.get() != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t10=\t*\t0\t0\tAAAAAAAAAA\t..........";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        BamAln aln("chr1");
        REQUIRE(aln.get_aln_internal() != nullptr);
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        CHECK(std::string(aln.get_name()) == "r1");
        CHECK(aln.get_ref_start() == 12344);
        CHECK(aln.get_ref_end() == 12353);
        CHECK(aln.get_query_start() == 0);
        CHECK(aln.get_query_len() == 10);
        CHECK(aln.get_cigar_len() == 1);
        CHECK(aln.get_seq().size() == 10);
        for (int i = 0; i < aln.get_seq().size(); i++) {
            CHECK(aln.get_seq()[i] == 'A');
        }
        CHECK(aln.get_mapq() == 60);
        CHECK(!aln.has_tag("SA"));
        CHECK(!aln.is_split());
        CHECK(!aln.is_reverse());
    }

    SUBCASE("split_single") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t0\tchr1\t12345\t60\t10=10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t....................\t"
                                 "SA:Z:chr1,23456,+,10S10M,40,0;";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        BamAln aln("chr1");
        REQUIRE(aln.get_aln_internal() != nullptr);
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        CHECK(std::string(aln.get_name()) == "r1");
        CHECK(aln.get_ref_start() == 12344);
        CHECK(aln.get_ref_end() == 12353);
        CHECK(aln.get_query_start() == 0);
        CHECK(aln.get_query_len() == 20);
        CHECK(aln.get_query_aln_len() == 10);
        CHECK(aln.get_cigar_len() == 2);
        CHECK(aln.get_seq().size() == 20);
        for (int i = 0; i < aln.get_seq().size(); i++) {
            CHECK(aln.get_seq()[i] == 'A');
        }
        CHECK(aln.get_mapq() == 60);
        CHECK(aln.has_tag("SA"));
        CHECK(aln.is_split());

        SplitAln split_aln = aln.make_split();
        CHECK(split_aln.chr_name == "chr1");
        CHECK(split_aln.ref_start == 12344);
        CHECK(split_aln.ref_end == 12353);
        CHECK(split_aln.query_start == 0);
        CHECK(split_aln.query_end == 10);
        CHECK(split_aln.mapq == 60);
        CHECK(!split_aln.is_reverse);

        std::vector<SplitAln> split_alns = aln.get_supplementary_alns();
        CHECK(split_alns.size() == 1);
        SplitAln& sa_aln = split_alns[0];
        CHECK(sa_aln.chr_name == "chr1");
        CHECK(sa_aln.ref_start == 23455);
        CHECK(sa_aln.ref_end == 23464);
        CHECK(sa_aln.query_start == 10);
        CHECK(sa_aln.query_end == 20);
        CHECK(sa_aln.mapq == 40);
        CHECK(!sa_aln.is_reverse);
    }

    SUBCASE("split_multi_no_gap") {
        REQUIRE(sam_hdr != nullptr);
        const char *sam_record = "r1\t16\tchr1\t12355\t60\t10S10=10S\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\t"
                                 "..............................\t"
                                 "SA:Z:chr1,12345,+,10M20S,20,0;chr1,12365,+,20S10M,10,0";
        kstring_t sam_entry = {0};
        kputs(sam_record, &sam_entry);
        BamAln aln("chr1");
        REQUIRE(aln.get_aln_internal() != nullptr);
        int ret = sam_parse1(&sam_entry, sam_hdr.get(), aln.get_aln_internal());
        CHECK(ret >= 0);
        CHECK(aln.get_ref_start() == 12354);
        CHECK(aln.get_ref_end() == 12363);
        CHECK(aln.get_query_start() == 10);
        CHECK(aln.get_query_len() == 30);
        CHECK(aln.get_query_aln_len() == 10);
        CHECK(aln.get_cigar_len() == 3);
        CHECK(aln.get_seq().size() == 30);
        for (int i = 0; i < aln.get_seq().size(); i++) {
            CHECK(aln.get_seq()[i] == 'A');
        }
        CHECK(aln.get_mapq() == 60);
        CHECK(aln.has_tag("SA"));
        CHECK(aln.is_split());
        CHECK(aln.is_reverse());

        std::vector<SplitAln> split_alns = aln.get_supplementary_alns();
        CHECK(split_alns.size() == 2);
        CHECK(split_alns[0].chr_name == "chr1");
        CHECK(split_alns[0].ref_start == 12344);
        CHECK(split_alns[0].ref_end == 12353);
        CHECK(split_alns[0].query_start == 0);
        CHECK(split_alns[0].query_end == 10);
        CHECK(split_alns[0].mapq == 20);
        CHECK(!split_alns[0].is_reverse);

        CHECK(split_alns[1].chr_name == "chr1");
        CHECK(split_alns[1].ref_start == 12364);
        CHECK(split_alns[1].ref_end == 12373);
        CHECK(split_alns[1].query_start == 20);
        CHECK(split_alns[1].query_end == 30);
        CHECK(split_alns[1].mapq == 10);
        CHECK(!split_alns[1].is_reverse);
    }
}