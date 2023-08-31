# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v0.10.0]
### Changed
Read iterator now returns copies of alignments for more Pythonic behaviour.

## [v0.9.5]
This version adds no user facing changes.

### Fixed
- Pinned pip in Makefile so CI tests will run.

## [v0.9.4]
### Fixed
- Python source distribution did not include libdeflate directory.

## [v0.9.3]
### Added
- Linux and macOS ARM conda builds.

## [v0.9.2]
### Changed
- Rebuild for conda with more specific htslib version pin.

## [v0.9.1]
### Added
- `--map_q` command line option to filter reads by minimum mapping quality.
### Fixed
- The default mapping quality was erroneously 1 not 0.

## [v0.9.0]
### Added
- "Other" modified base column to extended output. For example, when using
  `-m 5mC` to count 5-methylcytosine in reads, the "other" column will
  enumerate counts of other cytosine modifications present. When using the
  option `--combine` this column will contain zero (the counts being included
  in the modified base count).
- `--theshold` option to replace both `-a` and `-b`.
### Fixed
- In line with the "other" column, when not using the `--combine` option
  the potential presence of other modifications in the same family of the
  base requested is taking into account. This has an effect of distributing
  some previously erroneous "canonical/non-modified" calls to "other" and
  "filtered" counts.
### Removed
- The options `-a` and `-b` are deprecated. Instead use `--threshold`.

## [v0.8.0]
### Added
- `--combine` option to combine calls from all modified bases in a family.
  The previous behaviour was that the non-modified (canonical) count would
  have be incremented. For example when searching for 5mC modifications
  with `-m 5mC` and a 5hmC base was present, the read would contribute
  to the canonical count, not the modified count.

## [v0.7.0]
### Added
- `--pileup` option to output full raw base counts rather than BED methyl.
### Changed
- `-c` no longer synonym to `--cpg`.
- `?`-style MM subtags now handle correctly with "missing" entries being recorded
  as "no call" rather than implied canonical.
- extended output now includes a 15th column for "no call" bases.
### Fixed
- Links in README.

## [v0.6.3]
### Changed
- Bumped htslib version to version 1.16 for fixes to MM tag parsing/validation.
- Change conda build back to bioconda::htslib since we're using a released version.

## [v0.6.2]
### Fixed
- Off-by-one in pointless BED field.

## [v0.6.1]
### Fixed
- ChEBI codes not cast correctly in Python API.
### Added
- Support ambiguous modified bases as listed in HTS tags specification.

## [v0.6.0]
## Changed
- Sites with no coverage now report "nan" methylation frequency and score.
## Added
- Option `--aggregate` to pair information from two strands and output additional files.
- Support for ChEBI codes in C and Python pileup API.

## [v0.5.3]
## Added
- Python 3.9 and 3.10 wheel builds.

## [v0.5.2]
## Changed
- Use commit `e51f72f` of htslib for `?` and `.` parsing of Mm tag.

## [v0.5.1]
## Added
- `--max_depth` argument, and do not limit by default.
- `--chh` and `--chg` filter options.

## [v0.5.0]
### Changed
- Decouple file opening from read iteration.
- Move Python pileup function to method of ModBam class.

## [v0.4.6]
### Changed
- Reworked compilation to remove argparser from Python module.
### Fixed
- Memory leak in modbampy.

## [v0.4.5]
### Fixed
- Unmasking of reference sites (again).
### Added
- Option `--mask`/`-k` to respect reference soft-masking.

## [v0.4.4]
### Fixed
- Logic error in filtering CpG sites for masked bases.

## [v0.4.3]
### Changed
 - Update modbampy version to match C code.

## [v0.4.2]
### Changed
- Include soft-masked reference positions when performing CpG filtering.

## [v0.4.1]
### Fixed
- Inaccuracies in README.

## [v0.4.0]
### Fixed
- Python pileup access after addition of multi-BAM support.
### Added
- Additional properties to alignment objects in Python API.

## [v0.3.3]
### Changed
- conda build now uses htslib from bioconda.

## [v0.3.2]
### Changed
- Updated software build to use official version 1.14 htslib release.
- Reorganised and updated README.

## [v0.3.1]
### Changed
- Build conda package with explicit libdeflate version.

## [v0.3.0]
### Added
- Multiple BAM parsing to streamline interaction with data from Guppy/MinKNOW.
### Changed
- Reference file must now be given before list of BAM files on command-line.
### Fixed
- Removed some debugging text.

## [v0.2.2]
### Changed
- Updated README to note Python package.

## [v0.2.1]
### Fixed
- Incorrect processing of non-primary alignments in Python API.
### Added
- Add Python packages, available on PyPI.
### Changed
- Updated htslib to version from samtools/dev.


## [v0.2.0]
### Fixed
- Segmentation fault on exit caused by double free of faidx member.
### Added
- Python API to pileup and read-level parsing.

## [v0.1.1]
### Fixed
- Check input files are present and readable rather than segfaulting.
### Changed
- Clearer error messaging.


## [v0.1.0]

First release.
