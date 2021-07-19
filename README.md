Modified-base BAM to BedMethyl
------------------------------

Simple demonstration program of converting Modified-base BAM file to BedMethyl.

Uses `methylation` branch of htslib from [jkbonfield](https://github.com/jkbonfield).

    git clone --recursive <repository>
    make pileup
    ./pileup 4k_ecoli.bam ecoli1:1-5000000 ecoli.fasta 0.66 > counts.bed

The forth parameter is used to filter calls: calls with probability of being modified
in (1 - p, p) are disgarded, calls with P < 1 - p are deemed to be unmodified and
calls with P > p are modified.
