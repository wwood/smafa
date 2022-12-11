# A read aligner for aligned sequences

[![Travis](https://img.shields.io/travis/wwood/smafa.svg?style=flat-square)](https://travis-ci.org/wwood/smafa)

Smafa attempts to align or cluster pre-aligned biological sequences, handling
sequences which are all the same length. The main use case is through
[SingleM](https://github.com/wwood/singlem), although it can be used
independently without issue to search and cluster other pre-aligned sequences.

## Installation

### Installation through Cargo

Smafa can be installed in the usual way [Rust](http://rust-lang.org/) packages
are installed. After installing Rust, smafa can be installed using cargo:

```
cargo install smafa
```

## Usage

To run the aligner, first make a db with `smafa makedb` and then query that
database with `smafa query`. To see how to use these modes, use e.g. `smafa
makedb -h`.

## Help
If you have any questions or comments, please raise an issue on the GitHub
repository, or just email Ben Woodcroft.

## License
Smafa is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT. It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The code is available at [https://github.com/wwood/smafa](https://github.com/wwood/smafa).
