# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
