# `misasim`
Simulate a misassembly for a given fasta.

<table>
  <tr>
    <td>
      <figure float="left">
        <img align="middle" src="docs/imgs/chr16_MATERNAL_40000000-45000000_original.png" width="100%">
        <figcaption>Original <a href="https://github.com/logsdon-lab/NucFlag"><code>nucflag</code></a> plot</figcaption>
      </figure>
      <figure float="left">
        <img align="middle" src="docs/imgs/chr16_MATERNAL_40000000-45000000.png" width="100%">
        <figcaption>After introducing multiple misjoins</figcaption>
      </figure>
    </td>
  </tr>
</table>

```bash
# Delete 10 kbp of sequence and join the sequences (misjoin).
# Apply to every unique sequence in the assembly.
echo '[{"mtype": "misjoin", "number": 100, "length": 10000}]' > 10100.json
misasim multiple \
  --path 10100.json \
  --seed 10100 \
  --infile hg002v1.1.fasta.gz \
  --outfile hg002v1.1_misassembled.fa \
  --group-by "^(.*?)_(.*?)$" \
  --outbedfile regions.bed
```

```
chrom_og	chrom_og_st	chrom_og_end	name	score	strand	chrom_new_st	chrom_new_end	item_rgb
chr16_MATERNAL	41689710	41699710	misjoin	0	.	41259710	41259710	255,0,0
chr16_MATERNAL	42259214	42269214	misjoin	0	.	41819214	41819214	255,0,0
chr16_MATERNAL	43191922	43201922	misjoin	0	.	42741922	42741922	255,0,0
chr16_MATERNAL	43699445	43709445	misjoin	0	.	43239445	43239445	255,0,0
```

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
  multiple           Simulate multiple misassembly types from an input JSON file
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

### Input
* Fasta file
  * Sequences to introduce misassemblies into.
* BED3 file
  * Restrict introduced misassemblies to these regions
* (Optional) JSON file
  * Only applicable to `misasim multiple`.
  * List multiple misassemblies to introduce.
  * Duplicates will be provided different seeds. (+1 increment)

### Output
* Fasta file
  * Sequences with misassemblies.
* BED9 file
  * Regions in new assembly.

    |column|desc|
    |-|-|
    |chrom|Chromsome name|
    |chromStart|Chromsome start of misassembly/region in original assembly|
    |chromEnd|Chromsome end of misassembly/region in original assembly|
    |name|Misassembly type. Include good|
    |score|Always '0'|
    |strand|Always '.'|
    |thickStart|Chromsome start of misassembly/region in new assembly|
    |thickStart|Chromsome end of misassembly/region in new assembly|
    |itemRGB|Red (255,0,0) if a misassembly. White (0,0,0) otherwise|

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
-l 5000
```

#### Generate multiple misassemblies with the regions specified grouped by chromosome name.
```bash
./target/release/misasim multiple \
-i test/data/HG002.fa.gz \
-r test/data/region.bed \
-g "$(?<chr>.*?)_.*?$"
-p test/data/multiple.json
```

```json
[
    {
        "mtype": "misjoin",
        "number": 1,
        "length": 3
    },
    {
        "mtype": "inversion",
        "number": 1,
        "length": 3
    }
]
```
> `test/data/multiple.json`
