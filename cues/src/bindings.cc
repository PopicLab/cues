#include "cues/cues.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

namespace py = pybind11;
PYBIND11_MODULE(pycues, m) {
    py::class_<cues::RefSegment>(m, "RefSegment")
            .def(py::init<const std::string&, const cues::ref_pos_t, const cues::ref_pos_t>())
            .def_readwrite("chr", &cues::RefSegment::chr)
            .def_readwrite("start", &cues::RefSegment::start)
            .def_readwrite("end", &cues::RefSegment::end)
            .def("__str__", &cues::RefSegment::to_string);
    py::bind_vector<std::vector<cues::RefSegment>>(m, "RefSegmentList");

    py::enum_<cues::AlnStatistic>(m, "AlnStatistic")
            .value("kRD", cues::AlnStatistic::kRD)
            .value("kRDAll", cues::AlnStatistic::kRDAll)
            .value("kClipLeft", cues::AlnStatistic::kClipLeft)
            .value("kClipRight", cues::AlnStatistic::kClipRight)
            .value("kIns", cues::AlnStatistic::kIns)
            .value("kDiv", cues::AlnStatistic::kDiv);
    py::bind_map<cues::aln_stats_map_t>(m, "AlnStatsMap");
    py::bind_vector<std::vector<cues::aln_stats_map_t>>(m, "AlnStatsMapList");

    py::enum_<cues::JOrientation>(m, "JOrientation", py::arithmetic())
            .value("kLR", cues::JOrientation::kLR)
            .value("kRL", cues::JOrientation::kRL)
            .value("kLL", cues::JOrientation::kLL)
            .value("kRR", cues::JOrientation::kRR)
            .export_values();

    py::enum_<cues::JSource>(m, "JSource")
            .value("kCigar", cues::JSource::kCigar)
            .value("kSplitRead", cues::JSource::kSplitRead)
            .value("kSubFragPair", cues::JSource::kSubFragPair)
            .value("kClip", cues::JSource::kClip)
            .value("kIns", cues::JSource::kIns);

    py::class_<cues::Junction>(m, "Junction")
            .def(py::init<>())
            .def(py::init<std::string, std::string, cues::ref_pos_t, cues::ref_pos_t,
                         cues::JOrientation, cues::JSource, int>(),
                 py::arg("chr1"),
                 py::arg("chr2"),
                 py::arg("pos1"),
                 py::arg("pos2"),
                 py::arg("orientation"),
                 py::arg("evidence_type"),
                 py::arg("support"))
            .def_readwrite("chr1", &cues::Junction::chr1)
            .def_readwrite("chr2", &cues::Junction::chr2)
            .def_readwrite("pos1", &cues::Junction::pos1)
            .def_readwrite("pos2", &cues::Junction::pos2)
            .def_readwrite("source", &cues::Junction::evidence_type)
            .def_readwrite("orientation", &cues::Junction::orientation)
            .def_readwrite("support", &cues::Junction::support)
            .def("len", &cues::Junction::len, "Length of the spanned reference segment (only for intra-chromosomal junctions)")
            .def("__str__", &cues::Junction::to_string);
    py::bind_vector<std::vector<cues::Junction>>(m, "JunctionList");

    py::class_<cues::NovelSeq>(m, "NovelSeq")
            .def_readwrite("chr", &cues::NovelSeq::chr)
            .def_readwrite("ref_pos", &cues::NovelSeq::ref_pos)
            .def_readwrite("read_pos", &cues::NovelSeq::read_pos)
            .def_readwrite("len", &cues::NovelSeq::len)
            .def_readwrite("seq", &cues::NovelSeq::seq);

    py::class_<cues::AlnPatch>(m, "AlnPatch")
        .def_readwrite("aln_id", &cues::AlnPatch::aln_id)
        .def_readwrite("new_span", &cues::AlnPatch::new_span);

    py::class_<cues::BreakpointStore>(m, "BreakpointStore")
            .def_readwrite("junctions", &cues::BreakpointStore::junctions)
            .def_readwrite("insertions", &cues::BreakpointStore::insertions)
            .def_readwrite("aln_patches", &cues::BreakpointStore::aln_patches);

    ////////////////////
    // Configuration ///
    ////////////////////
    py::class_<cues::Cues::CuesConfig>(m, "CuesConfig")
            .def(py::init<>())
            .def_readwrite("ref", &cues::Cues::CuesConfig::ref, "Path to the reference genome FASTA")
            .def_readwrite("min_mapq", &cues::Cues::CuesConfig::min_mapq, "Minimum alignment MAPQ required for signal extraction")
            .def_readwrite("max_discordance", &cues::Cues::CuesConfig::max_discordance, "Maximum allowed alignment discordance")
            .def_readwrite("min_discord_len", &cues::Cues::CuesConfig::min_discord_len, "Minimum length of alignment discordance to record or to trigger signal recovery")
            .def_readwrite("max_discord_len", &cues::Cues::CuesConfig::max_discord_len, "Maximum size of novel junctions to record")
            .def_readwrite("min_overlap_thr", &cues::Cues::CuesConfig::min_overlap_thr, "Minimum overlap of split/fragment alignments required to assign kOrd21")
            .def_readwrite("min_accord_len", &cues::Cues::CuesConfig::min_accord_len, "Minimum length of alignment concordance to break discordant blocks")
            .def_readwrite("min_clipped_len", &cues::Cues::CuesConfig::min_clipped_len)
            .def_readwrite("max_segment_missing", &cues::Cues::CuesConfig::max_segment_missing)
            .def_readwrite("frag_len", &cues::Cues::CuesConfig::frag_len, "Length of fragments used for signal recovery")
            .def_readwrite("min_mapq_frag", &cues::Cues::CuesConfig::min_mapq_frag, "Minimum fragment alignment MAPQ required for junction detection")
            .def_readwrite("n_threads", &cues::Cues::CuesConfig::n_threads, "Number of threads for parallel processing across alignments/reads")
            .def_readwrite("aln_batch_size", &cues::Cues::CuesConfig::aln_batch_size, "Number of alignments to load at a time in each thread")
            .def_readwrite("max_alns_to_process", &cues::Cues::CuesConfig::max_alns_to_process, "Maximum number of alignments to process from BAM")
            .def_readwrite("max_seg_alns_to_fetch", &cues::Cues::CuesConfig::max_seg_alns_to_fetch, "Maximum number of alignments to fetch per segment")
            .def_readwrite("record_inter_chromosomal_junctions", &cues::Cues::CuesConfig::record_inter_chromosomal_junctions, "Collect inter-chromosomal junctions")
            .def_readwrite("record_insertions", &cues::Cues::CuesConfig::record_insertions, "Collect insertion records")
            .def_readwrite("record_aln_patches", &cues::Cues::CuesConfig::record_aln_patches, "Collect patched alignment records (enable for downstream stats reporting)")
            .def_readwrite("run_recovery_ins", &cues::Cues::CuesConfig::run_recovery_ins, "Run INS-based signal recovery")
            .def_readwrite("run_recovery_subfrags", &cues::Cues::CuesConfig::run_recovery_subfrags, "Run subfragment-based signal recovery")
            .def_readwrite("run_recovery_missing_segments", &cues::Cues::CuesConfig::run_recovery_missing_segments, "Run missing-segment signal recovery")
	    .def_readwrite("run_recovery_clipping", &cues::Cues::CuesConfig::run_recovery_clipping, "Run clipping-based signal recovery");

    py::class_<cues::Cues>(m, "Cues")
            .def(py::init<cues::Cues::CuesConfig&>())
            .def("collect_bkps_from_bam", &cues::Cues::collect_bkps_from_bam, py::call_guard<py::gil_scoped_release>(), "Collect breakpoint evidence from BAM")
            .def("collect_stats", &cues::Cues::collect_stats,
                     py::arg("bam_file_name"),
                     py::arg("segments"),
                     py::arg("stats_to_compute"),
                     py::arg("aln_patches") = std::vector<cues::AlnPatch>(),
                     "Compute alignment stats from BAM for a set of reference segments, given the pre-computed breakpoints");

}
