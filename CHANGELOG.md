# Changelog

## Unreleased

### Added

- Add entropy-balanced banding and pigeonhole prefilters for faster clustering and querying.
- Add multithreaded sequence processing.

### Changed

- Make entropy-balanced banding the default production prefilter.

## Version 0.8.0

### Added

- Add compact nucleotide encoding and support lowercase nucleotide symbols.
- Add database compatibility tests.

### Fixed

- Correct database version handling and improve invalid-nucleotide errors.

### Changed

- Update dependencies and improve command-line argument validation.

## Version 0.7.1

### Added

- Support degenerate nucleotide notation.

### Changed

- Update needletail to 0.5.

## Version 0.7.0

### Added

- Reintroduce clustering mode.

## Version 0.6.1

### Added

- Add the GPL license text.

### Changed

- Apply clippy fixes.

## Version 0.6.0

### Added

- Add count mode, top-hit queries, query result limits, CI, and release tooling.

### Changed

- Move sequence storage to serde/postcard and improve command-line help.

## Version 0.5.0

### Added

- Add amino-acid input, translated searches, and configurable k-mer length.

## Version 0.4.0

### Changed

- Increase divergence tolerance for fragment clustering and update dependencies.

## Version 0.3.0

### Fixed

- Fix duplicate output and incorrect fragment clustering.

## Version 0.2.0

### Added

- Add fragment-method clustering.

## Version 0.1.2

### Fixed

- Fix singleton detection and clustering output.

## Version 0.1.0

### Added

- Initial release.
