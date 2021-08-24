# COVID-Mixed-Lineage-Quantification
Quantify the abundance of different lineages in a mixed SARS-CoV-2 sequencing dataset

Minimap2 command:

```bash
minimap2 --MD -t 7 -a -x sr ~/Desktop/reference.fasta ~/Desktop/SEARCH-40821/*.fastq.gz | samtools view -bS - > SEARCH-40821.pangolin_ref.bam
```
