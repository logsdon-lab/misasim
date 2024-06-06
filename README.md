# `misasim`
Simulate a misassembly for a given fasta.


### Getting Started
Install [Rust](https://www.rust-lang.org/tools/install).

Compile `misasim`.
```bash
cargo build --release
```

### Usage
```
Usage: misasim <COMMAND>

Commands:
  misjoin            Simulate a misjoin in a sequence
  false-duplication  Simulate a falsely duplicated sequence
  gap                Simulate a gap in a sequence
  break              Simulate a break in a sequence
  help               Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
```

### Test
```bash
./target/release/misasim misjoin -i <(cat test/data/HG00171_chr9_haplotype2-0000142.fa) -b test.bed > test.fa
```

### TODO
* Add concurrent record reading with noodle async feature.
* Add gzip support.
* Write declarative macro to reduce DRY in cli.
