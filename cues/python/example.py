import argparse
import pycues

########## CLI ############
parser = argparse.ArgumentParser(description="A simple pycues example")
parser.add_argument("--ref", required=True, help="Path to the reference FASTA file")
parser.add_argument("--bam", required=True, help="Path to the BAM file")
args = parser.parse_args()
###########################

# initialize and populate the config
config = pycues.CuesConfig()
config.ref = args.ref
cues = pycues.Cues(config)

# extract breakpoints
breakpoints = cues.collect_bkps_from_bam(args.bam)
print("Number of novel adjacencies: ", len(breakpoints.junctions))
max_display = 10
for j in breakpoints.junctions[:max_display]:
    print("Junction: ", j.chr1, j.chr2, j.pos1, j.pos2, j.source, j.orientation, j.support)

# extract alignment stats at breakpoints
buffer = 1000
segments = []
segments = [pycues.RefSegment(j.chr1, j.pos1, j.pos1 + buffer) for j in breakpoints.junctions[:max_display]]

stats = {pycues.AlnStatistic.kRD,
         pycues.AlnStatistic.kRDAll,
         pycues.AlnStatistic.kClipLeft,
         pycues.AlnStatistic.kClipRight,
         pycues.AlnStatistic.kIns}
aln_stats = cues.collect_stats(args.bam, segments, stats)
for sid, stats_map in enumerate(aln_stats):
    for s in stats:
        print(sid, s, stats_map[s])
