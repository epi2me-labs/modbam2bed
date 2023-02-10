![Oxford Nanopore Technologies logo](https://github.com/epi2me-labs/modbam2bed/raw/master/images/ONT_logo_590x106.png)

Modified-base BAM to bedMethyl
------------------------------

A program to aggregate modified base counts stored in a
[modified-base BAM](https://samtools.github.io/hts-specs/SAMtags.pdf) (Section 2.1) file to 
a [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/) file.

A Python module is also available to obtain modified base information
from BAM files in a convenient form. It is envisaged that this will eventually
be replaced by an implementation in [pysam](https://pysam.readthedocs.io/en/latest/index.html).

### Installation

The program is available from our conda channel, so can be installed with:

    mamba create -n modbam2bed -c bioconda -c conda-forge -c epi2melabs modbam2bed

Packages are available for both Linux and MacOS.

Alternatively to install from the source code, clone the repository and then use make:

    git clone --recursive https://github.com/epi2me-labs/modbam2bed.git
    make modbam2bed
    ./modbam2bed

See the Makefile for more information. The code has been tested on MacOS (with
dependencies from brew) and on Ubuntu 18.04 and 20.04.

### Usage

The code requires aligned reads with the `Mm` and `Ml` tags (`MM` and `ML` also supported),
and the reference sequence used for alignment.

The below is a snapshot of the command-line interface; it may not be up-to-date, please
refer to the program `--help` option for the most accurate guidance.

```
Usage: modbam2bed [OPTION...] <reference.fasta> <reads.bam> [<reads.bam> ...]
modbam2bed -- summarise one or more BAM with modified base tags to bedMethyl. 

 General options:
      --aggregate            Output additional aggregated (across strand)
                             counts, requires --cpg or --chg.
      --combine              Create output with combined modified counts: i.e.
                             alternative modified bases within the same family
                             (same canonical base) are included.
  -c, --pileup               Output (full) raw base counts rather than BED
                             file.
  -e, --extended             Output extended bedMethyl including counts of
                             canonical, modified, and filtered bases (in that
                             order).
  -m, --mod_base=BASE        Modified base of interest, one of: 5mC, 5hmC, 5fC,
                             5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao. (Or modA,
                             modC, modG, modT, modU, modN for generic modified
                             base).
  -p, --prefix=PREFIX        Output file prefix. Only used when multiple output
                             filters are given.
  -r, --region=chr:start-end Genomic region to process.
  -t, --threads=THREADS      Number of threads for BAM processing.

 Base filtering options:
  -a, --canon_threshold=THRESHOLD
                             Deprecated. The option will be removed in a future
                             version. Please use --threshold.
  -b, --mod_threshold=THRESHOLD   Deprecated. The option will be removed in a
                             future version. Please use --threshold.
      --chg                  Output records filtered to CHG sites.
      --chh                  Output records filtered to CHH sites.
      --cpg                  Output records filtered to CpG sites.
  -f, --threshold=THRESHOLD  Bases with a call probability < THRESHOLD are
                             filtered from results (default 0.66).
  -k, --mask                 Respect soft-masking in reference file.

 Read filtering options:
  -d, --max_depth=DEPTH      Max. per-file depth; avoids excessive memory
                             usage.
  -g, --read_group=RG        Only process reads from given read group.
      --haplotype=VAL        Only process reads from a given haplotype.
                             Equivalent to --tag_name HP --tag_value VAL.
      --tag_name=TN          Only process reads with a given tag (see
                             --tag_value).
      --tag_value=VAL        Only process reads with a given tag value.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version

Mandatory or optional arguments to long options are also mandatory or optional
for any corresponding short options.
```

### Method and output format

Oxford Nanopore Technogies' sequencing chemistries and basecallers can detect
any number of modified bases. Compared to traditional methods which force a
false dichoctomy between say cytosine and 5-methylcytosine, this rich biology
needs to be remembered when interpreting modified base calls.

The htslib pileup API is used to create a matrix of per-strand base counts
including substitutions, modified bases and deletions. Inserted bases are not
counted. Bases of an abiguous nature (refered to as "filtered" below), as
defined by the filter threshold probabilities option `-b` are masked and used
(along with substitutions and deletions) in the definition of the "score"
(column 5) and "coverage" (column 10) entries of the bedMethyl file.

In the case of `?`-style `MM` subtags, where a lack of a recorded call should
not be taken as implying a canonical-base call, the "no call" count is incremented.
The "no call" count is used in the calculation of "coverage" and also the denominator
of "score".

In summary, a base is determined as being either "canonical", "modified", "filtered",
or "no call". The final output includes a modification frequency and score and
coverage information in order to assess the reliability of the frequency.

**Call filtering**

To determine the base present at a locus in a read, the query base in the
BAM record is examined along with the modified base information. A "canonical"
base probability is calculated as `1 - sum(P_mod)`, with `P_mod` being
the set of probabilities associated with all the modifications enumerated
in the BAM record. The base form with largest probability is taken as the
base present subject to the user-specified threshold. If the probability
is below the threshold the call is masked and contributes to the "filtered"
base count rather than the "canonical" or "modified" counts.

**Special Handling of alternative modified bases (`--combine` option)**

To intepret the case of multiple modifications being listed in
the BAM, `modbam2bed` can operate in two modes:

* *default*: alternative modified bases in the same family as the requested
  modification are counted separatedly as "other" --- neither in
  the "canonical" count of the "modified" count.
* `--combine`: alternative modified bases are lumped together into the 
  "modified" count and ultimately into a single modification frequency.

***A particular case where `--combine` is useful is when comparing to the result of bisulfite sequencing.***

**Output format**

> The description of the [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/)
> format on the ENCODE project website is rather loose. The definitions below are chosen pragmatically.

The table below describes precisely the entries in each column of the output BED
file. Columns seven to nine inclusive are included for compatibility with the BED
file specification, the values written are fixed and no meaning should be derived
from them. Columns 5, 10, and 11 are defined in terms of counts of observed
bases to agree with reasonable interpretations of the bedMethyl specifications:

 * N<sub>canon</sub> - canonical (unmodified) base count, (contigent on the use of `--combine`, see above.)
 * N<sub>mod</sub> - modified base count.
 * N<sub>filt</sub> - count of bases where read does not contain a substitution or deletion
   with respect to the reference, but the modification status is ambiguous: these bases
   were filtered from the calculation of the modification frequency.
 * N<sub>sub</sub> - count of reads with a substitution with respect to the reference.
 * N<sub>del</sub> - count of reads with a deletion with respect to the reference.
 * N<sub>no call</sub> - counts of reads with an absent modification call (but not a substitution or deletion).
 * N<sub>alt mod</sub> - counts of reads with and alternative modification call (but not a substitution or deletion).

Since these interpretations may differ from other tools an extended output is
available (enabled with the `-e` option) which includes three additional columns
with verbatim base counts.

| column | description                                                                                                                                                                                                                                                                  |
|--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1      | reference sequence name                                                                                                                                                                                                                                                      |
| 2      | 0-based start position                                                                                                                                                                                                                                                       |
| 3      | 0-based exclusive end position (invariably start + 1)                                                                                                                                                                                                                        |
| 4      | Abbreviated name of modified-base examined                                                                                                                                                                                                                                   |
| 5      | "Score" 1000 * (N<sub>mod</sub> + N<sub>canon</sub>) / (N<sub>mod</sub> + N<sub>canon</sub> + N<sub>no call</sub> + N<sub>alt mod</sub> + N<sub>filt</sub> + N<sub>sub</sub> + N<sub>del</sub>). The quantity reflects the extent to which the calculated modification frequency in Column 11 is confounded by the alternative calls. The denominator here is the total read coverage as given in Column 10. |
| 6      | Strand (of reference sequence). Forward "+", or reverse "-".                                                                                                                                                                                                                 |
| 7-9    | Ignore, included simply for compatibility.                                                                                                                                                                                                                                   |
| 10     | Read coverage at reference position including all canonical, modified, undecided (no calls and filtered), substitutions from reference, and deletions.  N<sub>mod</sub> + N<sub>canon</sub> + N<sub>no call</sub> + N<sub>alt mod</sub> + N<sub>filt</sub> + N<sub>sub</sub> + N<sub>del</sub>                                        |
| 11     | Percentage of modified bases, as a proportion of canonical and modified (excluding no calls, filtered, substitutions, and deletions).  100 \* N<sub>mod</sub> / (N<sub>mod</sub>  + N<sub>alt mod</sub> + N<sub>canon</sub>)                                                                                       |
| 12\*    | N<sub>canon</sub>                                                                                                                                                                                                                                                            |
| 13\*    | N<sub>mod</sub>                                                                                                                                                                                                                                                         |
| 14\*    | N<sub>filt</sub> those bases with a modification probability falling between given thresholds.                                                                                                                                                                           |
| 15\*    | N<sub>no call</sub> those bases for which the query base was the correct canonical base for the modified base being considered, but no call was made (see the definition of the `.` and `?` flags in the SAM tag specification).                                                                                                                                                                           |
| 16\*    | N<sub>alt mod</sub> those bases for which the query base was the correct canonical base for the modified base being considered, but and alternative modification was present.                                                                                                                                                                           |

\* Included in extended output only.


### Limitations

The code has not been developed extensively and currently has some limitations:

 * Support for motif filtering is limited to CpG, CHG, and CHH, sites. Without
   this filtering enabled all reference positions that are the canonical base
   (on forward or reverse strand) equivalent to the modified base under
   consideration are reported.
 * Insertion columns are completely ignored for simplicitly (and avoid
   any heuristics).
 * Second strand `MM` subtags (i.e. `MM:C-m` as compared with `MM:C+m`)
   are not supported. These are not typically used so shouldn't affect most users.
   If such a tag is detected and warning will be thrown and the tag ignored. These tags
   do come in to play for duplex basecalls.

### Python package

A Python package is available on [PyPI](https://pypi.org/project/modbampy/) which
contains basic functionality for parsing BAM files with modified-base information.
It is envisaged that this will eventually be replaced by an implementation in
[pysam](https://pysam.readthedocs.io/en/latest/index.html). As such the interface
is supplements but does not integrate or replace pysam.

The package can be installed with:

```
pip install modbampy
```

The package contains simply to modes of use. Firstly an interface to iterate
over reads in a BAM file and report modification sites:

```
from modbampy import ModBam
with ModBam(args.bam) as bam:
    for read in bam.reads(args.chrom, args.start, args.end):
        for pos_mod in read.mod_sites:
            print(*pos_mod)
```

Each line of the above reports the

* read_id,
* reference position,
* query (read) position,
* reference strand (+ or -),
* modification strand (0 or 1, as defined in the HTSlib tag specification. This is invariable 0),
* canonical base associated with modification,
* modified base,
* modified-base score (scaled to 0-255).

A second method is provided which mimics the couting procedure implemented in
`modbam2bed`:

```
from modbampy import ModBam
with ModBam(args.bam) as bam:
    positions, counts = bam.pileup(
        args.chrom, args.start, args.end
        low_threshold=0.33, high_threshold=0.66, mod_base="m")
```

The result is two [numpy](https://numpy.org/) arrays. The first indicates the reference
positions associated with the counts in the second array. Each row of the second array
(`counts` above) enumerates the observed counts of bases in the order:

    a c g t A C G T d D m M f F n N

where uppercase letters refer to bases on the forward strand, lowercase letters
relate to the reverse strand:

* A, C, G, T are the usual DNA bases,
* D indicates deletion counts,
* M modified base counts,
* F filtered counts - bases in reads with a modified-base record but which were filtered
  according to the thresholds provided.
* N no call base counts.

**Extras**

The read iterator API also contains a minimal set of functionality mirroring properties of 
alignments available from pysam. See the [code](https://github.com/epi2me-labs/modbam2bed/blob/master/modbampy/__init__.py)
for further details.

### Acknowledgements

We thank [jkbonfield](https://github.com/jkbonfield) for developing the modified base
functionality into the htslib pileup API, and [Jared Simpson](https://github.com/jts)
for testing and comparison to his independently developed code.

### Help

**Licence and Copyright**

© 2021- Oxford Nanopore Technologies Ltd.

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
