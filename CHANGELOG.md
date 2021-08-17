# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unrelease]
### Changed
- Updated htslib to version from samtools/dev.
### Fixed
- Incorrect processing of non-primary alignments in Python API.

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
