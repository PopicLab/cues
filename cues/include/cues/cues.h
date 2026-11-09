// Copyright (c) 2026 Victoria Popic
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <cstdint>
#include <format>
#include <limits>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifndef CUES_INCLUDE_CUES_CUES_H_
#define CUES_INCLUDE_CUES_CUES_H_

/**
 * Functionality for structural variation signal recovery
 * from read alignments, raw read sequences, or assemblies
 */
namespace cues {
    /////////////////
    // Basic types //
    /////////////////
    using ref_pos_t = int64_t; // reference position

    struct RefSegment { // interval in the reference genome
        std::string chr;
        ref_pos_t start; // 0-based, inclusive
        ref_pos_t end; // 0-based, exclusive
        RefSegment(const std::string& chr, const ref_pos_t start, const ref_pos_t end): chr(chr), start(start), end(end) {}
        inline ref_pos_t len() const { return abs(end - start); }
        inline std::string to_string() const { return std::format("{}:{}-{}", chr, start, end); };
    };

    /////////////////
    // SV signals  //
    /////////////////

    // Alignment statistics extracted from reads aligned to a reference locus/segment
    enum class AlnStatistic {
        kRD, // average read depth (alignments meeting the mapq threshold)
        kRDAll, // average read depth (all alignments)
        kClipLeft, // number of alignments clipped on the left
        kClipRight, // number of alignments clipped on the right
        kIns, // number of alignments containing an insertion
        kDiv, // average divergence (edits / reference span)
    };
    const std::unordered_map<AlnStatistic, std::string> STATS_TO_NAME = {
            {AlnStatistic::kRD, "RD"},
            {AlnStatistic::kRDAll, "RDAll"},
            {AlnStatistic::kClipLeft, "ClipLeft"},
            {AlnStatistic::kClipRight, "ClipRight"},
            {AlnStatistic::kIns, "Ins"},
            {AlnStatistic::kDiv, "Div"},
    };

    using aln_stats_value_t = double;
    using aln_stats_map_t = std::unordered_map<AlnStatistic, aln_stats_value_t>;

    // Rearrangement-induced novel junctions (defined by how pairs of consecutive genome segments map to the reference)
    enum class JOrientation { // segment orientation at the junction
        kLR,  // segments map to the same strand, in the order they appear in the read
        kRL,  // segments map to the same strand, out of order
        kLL, // segments map to two different strands, LL configuration of reads at each breakend
        kRR, // segments map to two different strands, RR configuration of reads at each breakend
    };

    const std::unordered_map<JOrientation, std::string> ORIENTATION_TO_NAME = {
            {JOrientation::kLR, "LR"},
            {JOrientation::kRL, "RL"},
            {JOrientation::kLL, "LL"},
            {JOrientation::kRR, "RR"},
    };

    enum class JSource { // source from which the junction was detected
        kCigar, // cigar (reported)
        kSplitRead, // split read alignment (reported)
        kSubFragPair,  // sub-fragment pair alignment (recovered)
        kClip, // clipped sequence alignment (recovered)
        kIns, // inserted sequence alignment (recovered)
    };

    const std::unordered_map<JSource, std::string> SOURCE_TO_NAME = {
            {JSource::kCigar, "Cigar"},
            {JSource::kSplitRead, "SplitRead"},
            {JSource::kSubFragPair, "SubFragPair"},
            {JSource::kClip, "Clip"},
            {JSource::kIns, "Ins"}
    };

    struct Junction {
        std::string chr1; // chromosome name of breakend 1
        std::string chr2; // chromosome name of breakend 2
        ref_pos_t pos1; // reference position of breakend 1
        ref_pos_t pos2; // reference position of breakend 2
        JOrientation orientation; // segment orientation at the junction
        JSource evidence_type; // source of evidence
        int support; // number of reads/alignments supporting this junction with matching evidence type, order, and orientation
        bool operator==(const Junction &) const = default;
        inline ref_pos_t len() const { return abs(pos1 - pos2); } // only for intra-chromosomal junctions
        inline std::string to_string() const {
            return std::format("{}@{}-{}@{}:{},{}", chr1, pos1, chr2, pos2,
                               ORIENTATION_TO_NAME.at(orientation),
                               SOURCE_TO_NAME.at(evidence_type));
        }
        struct Hash {
            size_t operator()(const Junction& j) const { return std::hash<std::string>()(j.to_string()); }
        };
    };
    using junctions_map_t = std::unordered_map<Junction, int, Junction::Hash>; // maps junctions to their occurrence

    struct NovelSeq { // novel insertion
        std::string chr;
        ref_pos_t ref_pos;
        ref_pos_t read_pos;
        ref_pos_t len;
        std::string seq;
    };

	struct AlnPatch {
		std::string aln_id;
		RefSegment new_span;
	};

	struct BreakpointStore {
        std::vector<Junction> junctions;
        std::vector<NovelSeq> insertions;
		std::vector<AlnPatch> aln_patches;
    };

    class Cues {
    public:
        // Parameters controlling SV signal extraction and execution
        struct CuesConfig {
            std::string ref; // path to the reference genome FASTA
            int min_mapq = 10; // minimum alignment MAPQ required for signal extraction
            int min_discord_len = 40; // minimum length of alignment discordance to record or to trigger signal recovery
            int max_discord_len = 1000000000; // maximum size of novel junctions to record
            int min_overlap_thr = 50; // minimum overlap of split/fragment alignments required to assign ORD21
            int min_accord_len = 25; // minimum length of alignment concordance to break discordant blocks
            int min_clipped_len = 10;
            int max_segment_missing = 15;
            double max_discordance = 1.0;
            int frag_len = 50; // length of fragments used for signal recovery
            int min_mapq_frag = 10; // minimum fragment alignment MAPQ required for junction detection
            int n_threads = 1; // number of threads for parallel processing across alignments/reads
            int aln_batch_size = 500000; // number of alignments/reads to load at a time in each thread
            int max_alns_to_process = std::numeric_limits<int>::max(); // maximum number of alignments to process from BAM
            int max_seg_alns_to_fetch = 100; // maximum number of alignments to fetch per segment
            bool record_inter_chromosomal_junctions = false;
            bool record_insertions = false;
            bool record_aln_patches = false;
            bool run_recovery_ins = false;
            bool run_recovery_clipping = false;
            bool run_recovery_subfrags = false;
            bool run_recovery_missing_segments = false;
        };
        explicit Cues(const CuesConfig& config): config_(config) {}
        inline const CuesConfig& get_config() { return config_; }
        // Extract breakpoints (novel adjacencies and insertions) from BAM
        BreakpointStore collect_bkps_from_bam(const std::string& bam_file_name);
        // Compute alignment stats from BAM for a set of reference segments, given the pre-computed breakpoints
        std::vector<aln_stats_map_t> collect_stats(const std::string& bam_file_name,
                                                   const std::vector<RefSegment>& segments,
                                                   const std::unordered_set<AlnStatistic>& stats_to_compute,
												   const std::vector<AlnPatch>& aln_patches = {});
    private:
        CuesConfig config_;
    };
}

#endif //CUES_INCLUDE_CUES_CUES_H_
