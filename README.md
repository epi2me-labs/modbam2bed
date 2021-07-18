Modified-base BAM to BedMethyl
------------------------------

Simple demonstration program of converting Modified-base BAM file to BedMethyl.

Uses `methylation` branch of htslib from (jkbonfield)[https://github.com/jkbonfield].

    git clone --recursive <repository>
    make pileup
    ./pileup 4k_ecoli.bam ecoli1:1-5000000 ecoli.fasta > counts.bed
