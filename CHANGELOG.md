# Changelog

All notable changes to codonyat are documented in this file.

## [1.2.0] - 2026-10-13

### Added
- **Multi-record reference FASTA support**: Load references with multiple sequence records (e.g., segmented viruses like Influenza A). Proteins are addressed by (record_id, protein_name) tuple.
- **Multi-interval proteins**: Support for complex protein structures:
  - Programmed ribosomal frameshift (e.g., SARS-CoV-2 ORF1ab with -1 PRF at base 13468, where base 13468 is read twice)
  - Splicing (discontinuous coding sequences)
  - Circular origin-crossing (e.g., HBV polymerase spanning the origin)
- **TSV annotation format**: Optional sidecar annotation file with format `record, protein, description, intervals, strand`
  - Intervals notation: 1-based inclusive coordinates, semicolon-separated for multi-interval proteins
  - Example: `266-13468;13468-21555` for SARS-CoV-2 ORF1ab frameshift
- **RECORD column in output**: For multi-record references, TSV output includes a `RECORD` column after `REFERENCE`
  - Single-record outputs remain byte-identical to v1.1.0 (no RECORD column)
- **Minus-strand gene validation**: Clear error message when minus-strand genes are encountered (not supported)
- **CLI flag**: `--annotation` for specifying TSV annotation files
- **Comprehensive test suite**: 37 tests covering all new features and backward compatibility

### Changed
- `FullReference` class now supports multiple records and optional annotation files
- `Protein` model now uses `Interval` objects and supports multiple intervals per protein
- `calculate_variant_frequencies` method updated to handle multi-interval proteins
- `write_csv` method conditionally includes RECORD column based on number of records

### Fixed
- Import formatting for better code organization
- Validation messages updated for multi-record references

### Backward Compatibility
- Single-record, single-interval workflows remain 100% byte-identical to v1.1.0
- Existing FASTA header annotation format continues to work unchanged
- All existing CLI flags and behavior preserved
- 26 original tests continue to pass unchanged

### Documentation
- Updated README with multi-interval protein examples
- Documented intervals notation: 1-based inclusive, semicolon-separated
- Added TSV annotation format specification
- Clarified RECORD column behavior
- Added "Research use only" wording throughout

## [1.0.3] - Previous Release

Previous stable version. See git history for details.
