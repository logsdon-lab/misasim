#!/bin/bash

set -euo pipefail


original_bam="/project/logsdon_shared/projects/HGSVC3/CenMAP_verkko_run2/results/nucflag/HG00171_hifi.bam"
reads_list="HG00171_chr9_haplotype2-0000142_reads.list"
reads_dir="/project/logsdon_shared/data/hgsvc3/hifi_reads_batch123/HG00171"
region_name="HG00171_chr9_haplotype2-0000142"

samtools view $original_bam "$region_name" | \
    cut -f 1 | \
    sort | \
    uniq > $reads_list

seqtk subseq $reads_dir/m54329U_220205_003428.hifi_reads.fastq.gz $reads_list | bgzip > m54329U_220205_003428.hifi_reads.fastq.gz
seqtk subseq $reads_dir/m54329U_220212_223600.hifi_reads.fastq.gz $reads_list | bgzip > m54329U_220212_223600.hifi_reads.fastq.gz
seqtk subseq $reads_dir/m54329U_220214_093304.hifi_reads.fastq.gz $reads_list | bgzip > m54329U_220214_093304.hifi_reads.fastq.gz