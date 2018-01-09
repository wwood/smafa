# A read aligner for aligned sequences

[![Travis](https://img.shields.io/travis/wwood/smafa.svg?style=flat-square)](https://travis-ci.org/wwood/smafa)

Smafa attempts to align or cluster pre-aligned biological sequences, handling
sequences which are all the same length. The main use case is through
[SingleM](https://github.com/wwood/singlem), although it can be used without
independently without issue.

## Installation

To run it, you'll need [Rust](http://rust-lang.org/) 1.20+. Then smafa can be
installed using cargo:

```
cargo install smafa
```

## Usage

To run the aligner, first make a db with `smafa makedb` and then query that
database with `smafa query`. To see how to use these modes, use e.g. `smafa
makedb -h`.

To run the in clustering mode, use `smafa cluster`. Note that the clustering
mode implements a greedy algorithm, where sequences encountered earlier in the
input file are taken as cluster representatives, unless they are sufficiently
similar to, i.e. cluster with, a previously encountered sequence.

## Help
If you have any questions or comments, send a message to the 
[SupportM mailing list](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/supportm)
or raise a [GitHib issue](https://github.com/wwood/smafa/issues).

## Development notes
To run benchmarks, rust nightly is required. Then run
```
rustup run nightly cargo bench --features unstable
```

## Citation
Woodcroft, B. Community diversity in metagenomes: one, many and thousands.
Winter School in Mathematical and Computational Biology.
http://bioinformatics.org.au/ws16/speaker-items/ben-woodcroft/#tab-f2e1404e449518dbab9

## License
Smafa is written by [BenWoodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft)
(@wwood) at the [Australian Centre for Ecogenomics (UQ)](http://ecogenomic.org/)
and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
