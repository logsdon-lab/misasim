```bash
conda activate general
python workflow/scripts/concensus.py         -i results/misasim/repair/HG00171_chr9_haplotype2-0000142_misjoin_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd_format.paf -r results/misasim/HG00171_chr9_haplotype2-0000142_misjoin.fa -q results/misasim/HG00171_chr9_haplotype2-0000142_false-duplication.fa         -b results/misasim/concensus/HG00171_chr9_haplotype2-0000142_misjoin_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd.bed -o results/misasim/concensus/HG00171_chr9_haplotype2-0000142_misjoin_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd.fa         --input_ref_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_misjoin_fwd.bed         --input_query_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_false-duplication_fwd.bed         --merge
```

```bash
python workflow/scripts/concensus.py         -i results/misasim/repair/HG00171_chr9_haplotype2-0000142_false-duplication_fwd_HG00171_chr9_haplotype2-0000142_misjoin_fwd_format.paf -r results/misasim/HG00171_chr9_haplotype2-0000142_false-duplication.fa -q results/misasim/HG00171_ch
r9_haplotype2-0000142_misjoin.fa         -b results/misasim/concensus/HG00171_chr9_haplotype2-0000142_false-duplication_fwd_HG00171_chr9_haplotype2-0000142_misjoin_fwd.bed -o results/misasim/concensus/HG00171_chr9_haplotype2-0000142_false-duplication_fwd_HG00171_chr9_haplotype2-0000142_mi
sjoin_fwd.fa         --input_ref_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_false-duplication_fwd.bed         --input_query_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_misjoin_fwd.bed         --merge
```

```bash
python workflow/scripts/concensus.py         -i results/misasim/repair/HG00171_chr9_haplotype2-0000142_gap_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd_format.paf -r results/misasim/HG00171_chr9_haplotype2-0000142_gap.fa -q results/misasim/HG00171_chr9_haplotype2-0000
142_false-duplication.fa         -b results/misasim/concensus/HG00171_chr9_haplotype2-0000142_gap_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd.bed -o results/misasim/concensus/HG00171_chr9_haplotype2-0000142_gap_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd.fa
     --input_ref_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_gap_fwd.bed         --input_query_misasm_bed results/misasim/concensus/HG00171_chr9_haplotype2-0000142_false-duplication_fwd.bed         --merge 2> logs/generate_concensus_HG00171_chr9_haplotype2-0000142
_gap_fwd_HG00171_chr9_haplotype2-0000142_false-duplication_fwd.log
```
