# `misasim`
Simulate a misassembly for a given fasta.

<table>
  <tr>
    <td>
      <figure float="left">
        <img align="middle" src="docs/imgs/HG00171_chr9_haplotype2-0000142:88082882-90122460_original.png" width="100%">
        <figcaption>Original</figcaption>
      </figure>
      <figure float="left">
        <img align="middle" src="docs/imgs/HG00171_chr9_haplotype2-0000142:88082882-90122460_misjoin.png" width="100%">
        <figcaption>After misjoin of 5000 bp</figcaption>
      </figure>
    </td>
  </tr>
</table>

### Getting Started
Install [Rust](https://www.rust-lang.org/tools/install).

Compile `misasim`.
```bash
cargo build --release
```

### Usage
```
Usage: misasim [OPTIONS] <COMMAND>

Commands:
  misjoin            Simulate a misjoin in a sequence
  false-duplication  Simulate a falsely duplicated sequence
  gap                Simulate a gap in a sequence
  break              Simulate a break in a sequence
  inversion          Simulate an inversion in a sequence
  help               Print this message or the help of the given subcommand(s)

Options:
  -i, --infile <INFILE>          Input sequence file. Uncompressed or bgzipped
  -r, --inbedfile <INBEDFILE>    Input bed file. Each region should map to a sequence from infile
  -o, --outfile <OUTFILE>        Output sequence file
  -b, --outbedfile <OUTBEDFILE>  Output BED file with misassemblies
  -s, --seed <SEED>              Seed to use for the random number generator
      --randomize-length         Randomize length
  -g, --group-by <GROUP_BY>      Group by regex pattern. ex. "^.*?_(?<hap>.*?)$" with group by haplotype
  -h, --help                     Print help
```

### Examples:

#### Generate a misjoin at a random position with a random length.
```bash
./target/release/misasim misjoin \
-i test/data/HG002_chr10_cens.fa.gz
```

#### Generate 12 misjoins at random positions with a random length.
```bash
./target/release/misasim misjoin \
-i test/data/HG002_chr10_cens.fa.gz \
-n 12
```

#### Generate a false-duplication at a random position with a length of 5000 bp.
```bash
./target/release/misasim false-duplication \
-i test/data/HG002_chr10_cens.fa.gz \
-l 5000
```

#### Generate a false-duplication at a random position with a length of 5000 bp duplicated at most four times.
```bash
./target/release/misasim false-duplication \
-i test/data/HG002_chr10_cens.fa.gz \
-l 5000 \
--max-duplications 4
```

#### Generate a gap at a random position within the regions specified with a length of 5000 bp.
```bash
./target/release/misasim misjoin \
-i test/data/HG002.fa.gz \
-r test/data/region.bed \
-l 5000
```

#### Generate a misjoin at a random position within the regions specified with a length of 5000 bp grouped by chromosome name.
```bash
# Either chr10_MATERNAL or chr10_PATERNAL would get a misjoin.
./target/release/misasim misjoin \
-i test/data/HG002.fa.gz \
-r test/data/region.bed \
-g "$(?<chr>.*?)_.*?$" # "$(.*?)_.*?$" would also work.
```
