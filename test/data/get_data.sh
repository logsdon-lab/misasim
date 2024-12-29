#!/bin/bash

set -euo pipefail

url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz"
fa="HG002.fa.gz"
bed="region.bed"
outfa="HG002_chr10_cens.fa.gz"

wget -O "${fa}" "${url}"
samtools faidx "${fa}"
seqtk subseq "${fa}" "${bed}" | bgzip > "${outfa}"
samtools faidx ${outfa}
