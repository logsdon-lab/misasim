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

```bash
./target/release/misasim misjoin -i test/data/HG00171_chr9_haplotype2-0000142.fa -b test.bed > test.fa
```

```bash
cat test/data/HG00171_chr9_haplotype2-0000142.fa | ./target/release/misasim misjoin -i - -b test.bed > test.fa
```

### Test
To test [`NucFlag`](https://github.com/logsdon-lab/NucFlag) on the misassembled contigs.
* This requires the original assembly and pacbio hifi reads used to produce it.
```
test/data/reads
├── m54329U_220205_003428.hifi_reads.fastq
├── m54329U_220212_223600.hifi_reads.fastq
└── m54329U_220214_093304.hifi_reads.fastq
test/data/asm/
├── HG00171_asm.fa
└── HG00171_asm.fa.fai
```

To run it on your cluster, `conda` and `Snakemake==7.32.4` are required.
```bash
queue_name="epistasis_normal"
snakemake -j 100 \
--cluster "bsub -q ${queue_name} -n {threads} -M {resources.mem}GB -R 'rusage[mem={resources.mem}GB]' -o /dev/null" \
-p \
-s test/workflow/Snakefile \
--use-conda \
-d test/
```

### TODO
* Add concurrent record reading with noodle async feature.
* Add gzip support.
* Write declarative macro to reduce DRY in cli.
