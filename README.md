Modified-base BAM to bedMethyl
------------------------------

Simple demonstration program of aggregating modified base counts stored in a
[modified-base BAM](https://circle-production-customer-artifacts.s3.amazonaws.com/picard/548f4caff7d0cea1406e35e3/60e6b57a449db620ab0ea9b7-0-build/artifacts/pdfs/SAMtags.pdf?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20210725T133456Z&X-Amz-SignedHeaders=host&X-Amz-Expires=60&X-Amz-Credential=AKIAJR3Q6CR467H7Z55A%2F20210725%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=0590263dd0113d97f83c533f1e3633cdba6a9652a8a9f457d6a0eb0d188037ef) (Section 2.1) file to 
a [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file.

The code uses the `methylation` branch of htslib from
[jkbonfield](https://github.com/jkbonfield).

### Installation

Clone the repository and then use make:

    git clone --recursive <repository>
    make modbam2bed
    ./modbam2bed

See the Makefile for more information. The code has been tested on MacOS (with
dependencies from brew) and on Ubuntu 18.04.

### Usage

The code requires aligned reads with the `Mm` and `Ml` tags, and the reference
sequence used for alignment.

```
Usage: modbam2bed [OPTION...]  <reads.bam> <reference.fasta> >
modbambed -- summarise a BAM with modified base tags to bedMethyl.

  -a, --canon_threshold=CANON_THRESHOLD
                             Bases with mod. probability < CANON_THRESHOLD are
                             counted as canonical.
  -b, --mod_threshold=MOD_THRESHOLD
                             Bases with mod. probability > MOD_THRESHOLD are
                             counted as modified.
  -e, --extended             Output extended bedMethyl including counts of
                             canonical, modified, and filtered bases (in that
                             order).
  -g, --read_group=READ_GROUP   Only process reads from given read group.
  -t, --threads=THREADS      Number of threads for BAM processing.
  -r, --region=chr:start-end Genomic region to process.
  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.

Positions absent from the methylation tags are assumed to be canonical.
Positions with modified probability between upper and lower thresholds are
removed from the counting process. Column 5 ("score") of the output is
calculated as the proportion of bases called as the canonical or modified
reference dase with respect to the number of spanning reads, scaled to a
maximum of 1000. Column 11 is the percentage of reference base calls identified
as being modified.

Report bugs to chris.wright@nanoporetech.com.
```

### Method

The htlib pileup API is used to create a matrix of per-strand base counts
including modified bases and deletions. Inserted bases are not counted. Bases
of an abiguous nature, as defined by the two threshold probabilities are
masked and used (along with substitutions and deletions) in the definition
of the "score" column of the bedMethyl file.

> The description of the bedMethyl format on the ENCODE project website is rather
> loose. The definitions below are chosen pragmatically.

The output bedMethyl file contains per-strand modified base proportions in column 11:

    100 * N_mod / (N_mod + N_canon)

where only canonical bases equivalent to the modified base are counted
(substitutions and deletions are ignored). Similarly column 5 contains the modified
count scaled by the total depth on the strand:

    1000 * (N_mod + N_canon) / (N_mod + N_canon + N_filt + N_sub + N_del)

Note the denominator here is the total sequencing depth at the position
(including deleted and ambiguous bases); the quantity therefore reflects the
extent to which the decision of modified or canonical base is confounded by
alternative calls.

#### Extended Output

As the bedMethyl definition is rather loose, the program allows for an "extended"
output which includes three additional columns reporting counts of relevant bases:

 * column 12: count of canonical base,
 * column 13: count of modified base,
 * column 14: count of filtered bases (those with a modification probability
   falling between the `--low_threshold` and `--high_threshold` parameters).

### Limitations

The code has not been developed extensively and currently has some limitations:

 * Only 5mC is currently supported as a modified base.
 * All reference positions that are C (on forward or reverse strand) are
   reported; there's no support for filtering by biologically interesting
   motifs.
 * No option to combine per-strand counts into a total count (how to do this
   generically depends on motif).
 * Insertion columns are completely ignored for simplicitely (and avoid
   any heuristics).

### Acknowledgements

We thank [jkbonfield](https://github.com/jkbonfield) for developing the modified base
functionality into the htslib pileup API, and [Jared Simpson](https://github.com/jts)
for testing and comparison to his independently developed code.

### Help

**Licence and Copyright**

© 2021 Oxford Nanopore Technologies Ltd.

`modbam2bed` is distributed under the terms of the Mozilla Public License 2.0.

**Research Release**

Research releases are provided as technology demonstrators to provide early
access to features or stimulate Community development of tools. Support for
this software will be minimal and is only provided directly by the developers.
Feature requests, improvements, and discussions are welcome and can be
implemented by forking and pull requests. However much as we would
like to rectify every issue and piece of feedback users may have, the
developers may have limited resource for support of this software. Research
releases may be unstable and subject to rapid iteration by Oxford Nanopore
Technologies.
