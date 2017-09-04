# A read aligner for aligned sequences

smafa attempts to align sequences to a database of aligned sequences. Currently,
it only handles query and db sequences which are all the same length, in concert
with [SingleM](https://github.com/wwood/singlem). To run it, you'll need
[Rust](http://rust-lang.org/) 1.20+. Then download and run smafa:

```
$ git clone https://github.com/wwood/smafa
$ cd smafa
$ cargo run --release -- -h

Ben J. Woodcroft <benjwoodcroft near gmail.com>
Read aligner for small pre-aligned sequences

USAGE:
    smafa [FLAGS] [OPTIONS] <DB_FASTA> <QUERY_FASTA>

FLAGS:
    -h, --help       Prints help information
    -q, --quiet      Unless there is an error, do not print logging information
    -V, --version    Prints version information
    -v, --verbose    Print extra debug logging information

OPTIONS:
    -d, --divergence <INTEGER>    Maximum number of mismatches in reported hits [default: 5]

ARGS:
    <DB_FASTA>       Subject sequences to search against
    <QUERY_FASTA>    Query sequences to search with
```

# License
Copyright Ben Woodcroft. Licensed under GPL3+. See LICENSE.txt in this
repository.
