# cues

`cues` is a C++ library (with Python bindings) for discovering structural variant (SV) signals
in long-read alignments. Given an input BAM file, `cues` extracts: (1) novel sequence adjacencies 
(junctions): i.e., pairs of reference positions that are distant in the reference but joined 
in the reads; and (2) alignment statistics over reference segments (e.g., read depth and clipping). 
`cues` reports both the junctions directly encoded in the input alignments 
(e.g., large CIGAR deletions and split reads) and junctions recovered from alignment artifacts
(using local realignment of read fragments).
Each reported junction is annotated with its orientation, evidence source (e.g., CIGAR),
and number of supporting alignments.

`cues` is based on the SV signal discovery approach developed in the
[Cue2](https://github.com/PopicLab/cue2) SV calling framework.

## Install

```bash
conda env create -f install/conda_env.yml
conda activate cues
./build.sh                 # C++ library and example (build/)
pip install ./cues         # Python module (pycues)
```

## Usage

Examples of how to use the library are provided for C++ in `cues/src/main.cc`
and for Python in `cues/python/example.py`.

```bash
./build/bin/cues_example <ref.fa> <aln.bam>
python cues/python/example.py --ref <ref.fa> --bam <aln.bam>
```

## Input requirements

`cues` requires two main inputs:

* BAM: sorted and indexed; supplementary alignments must be soft-clipped (no hard clipping);
  for best results, CIGAR strings should use `=`/`X` operations.
* Reference FASTA (indexed) used for alignment.

## Authors

Victoria Popic (vpopic@broadinstitute.org)
