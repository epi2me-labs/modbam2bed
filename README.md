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

```
Usage: modbam2bed [OPTION...] <reference.fasta> <reads.bam> [<reads.bam> ...]
modbam2bed -- summarise one or more BAM with modified base tags to bedMethyl.

 General options:
  -e, --extended             Output extended bedMethyl including counts of
                             canonical, modified, and filtered bases (in that
                             order).
  -m, --mod_base=BASE        Modified base of interest, one of: 5mC, 5hmC, 5fC,
                             5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao.
  -p, --prefix=PREFIX        Output file prefix. Only used when multiple output
                             filters are given.
  -r, --region=chr:start-end Genomic region to process.
  -t, --threads=THREADS      Number of threads for BAM processing.

 Base filtering options:
  -a, --canon_threshold=THRESHOLD
                             Bases with mod. probability < THRESHOLD are
                             counted as canonical (default 0.33).
      --aggregated           Output additional aggregated (across strand)
                             counts, requires --cpg or --chg.
  -b, --mod_threshold=THRESHOLD   Bases with mod. probability > THRESHOLD are
                             counted as modified (default 0.66).
      --cpg                  Output records filtered to CpG sites.
      --chg                  Output records filtered to CHG sites.
      --chh                  Output records filtered to CHH sites.
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
```

### Method and output format

The htslib pileup API is used to create a matrix of per-strand base counts
including modified bases and deletions. Inserted bases are not counted. Bases
of an abiguous nature, as defined by the two threshold probabilities are
masked and used (along with substitutions and deletions) in the definition
of the "score" (column 5) and "coverage" (column 10) entries of the bedMethyl file.

> The description of the [bedMethyl](https://www.encodeproject.org/data-standards/wgbs/)
> format on the ENCODE project website is rather loose. The definitions below are chosen pragmatically.

The table below describes precisely the entries in each column of the output BED
file. Columns seven to nine inclusive are included for compatibility with the BED
file specification, the values written are fixed and no meaning should be derived
from them. Columns 5, 10, and 11 are defined in terms of counts of observed
bases to agree with reasonable interpretations of the bedMethyl specifications:

 * N<sub>canon</sub> - canonical (unmodified) base count.
 * N<sub>mod</sub> - modified base count.
 * N<sub>filt</sub> - count of bases where read does not contain a substitution or deletion
   with respect to the reference, but the modification status is ambiguous: these bases
   were filtered from the calculation of the modification frequency.
 * N<sub>sub</sub> - count of reads with a substitution with respect to the reference.
 * N<sub>del</sub> - count of reads with a deletion with respect to the reference.

Since these interpretations may differ from other tools an extended output is
available (enabled with the `-e` option) which includes three additional columns
with verbatim base counts.

| column | description                                                                                                                                                                                                                                                                  |
|--------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1      | reference sequence name                                                                                                                                                                                                                                                      |
| 2      | 0-based start position                                                                                                                                                                                                                                                       |
| 3      | 0-based exclusive end position (invariably start + 1)                                                                                                                                                                                                                        |
| 4      | Abbreviated name of modified-base examined                                                                                                                                                                                                                                   |
| 5      | "Score" 1000 * (N<sub>mod</sub> + N<sub>canon</sub>) / (N<sub>mod</sub> + N<sub>canon</sub> + N<sub>filt</sub> + N<sub>sub</sub> + N<sub>del</sub>). The quantity reflects the extent to which the calculated modification frequency in Column 11 is confounded by the alternative calls. The denominator here is the total read coverage as given in Column 10. |
| 6      | Strand (of reference sequence). Forward "+", or reverse "-".                                                                                                                                                                                                                 |
| 7-9    | Ignore, included simply for compatibility.                                                                                                                                                                                                                                   |
| 10     | Read coverage at reference position including all canonical, modified, undecided (filtered), substitutions from reference, and deletions.  N<sub>mod</sub> + N<sub>canon</sub> + N<sub>filt</sub> + N<sub>sub</sub> + N<sub>del</sub>                                        |
| 11     | Percentage of modified bases, as a proportion of canonical and modified (excluding filtered, substitutions, and deletions).  100 \* N<sub>mod</sub> / (N<sub>mod</sub> + N<sub>canon</sub>)                                                                                       |
| 12\*    | N<sub>canon</sub>                                                                                                                                                                                                                                                            |
| 13\*    | N<sub>mod</sub>                                                                                                                                                                                                                                                         |
| 14\*    | N<sub>filt</sub> those bases with a modification probability falling between given thresholds.                                                                                                                                                                           |

\* Included in extended output only.


### Limitations

The code has not been developed extensively and currently has some limitations:

 * Support for motif filtering is limited to CpG, CHG, and CHH, sites. Without
   this filtering enabled all reference positions that are the canonical base
   (on forward or reverse strand) equivalent to the modified base under
   consideration are reported.
 * Insertion columns are completely ignored for simplicitly (and avoid
   any heuristics).
 * Results for opposite strand `MM` tags (i.e. `MM:C-m` as compared with `MM:C+m`)
   is not well tested. These are not typically used so shouldn't affect most users.
   They do come in to play for duplex basecalls.

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

    a c g t A C G T d D m M f F

where uppercase letters refer to bases on the forward strand, lowercase letters
relate to the reverse strand:

* A, C, G, T are the usual DNA bases,
* D indicates deletion counts,
* M modified base counts,
* F filtered counts - bases in reads with a modified-base record but which were filtered
  according to the thresholds provided.

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
