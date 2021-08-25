# CovidQuant
Quantify the abundance of different lineages in a mixed SARS-CoV-2 sequencing dataset

## Convert Trimmed BAM to Trimmed FASTQ

```bash
bedtools bamtofastq -i <TRIMMED_BAM> -fq /dev/stdout | pigz -9 > <TRIMMED_FASTQ_GZ>
```

Batch:

```bash
for f in *.bam; do bedtools bamtofastq -i "$f" -fq /dev/stdout | pigz -9 > "$(echo $f | cut -d'_' -f1).trimmed.fastq.gz" ; done
```

## Remap Using Minimap2

Need to remap using [Pangolin reference genome](https://github.com/cov-lineages/pangolin/blob/master/pangolin/data/reference.fasta).

```bash
minimap2 -a -x sr <REFERENCE_FASTA> <TRIMMED_FASTQ_GZ> 2> <MINIMAP2_LOG> | samtools view -bS - > <TRIMMED_REALIGNED_BAM>
```

Batch:

```bash
for f in *.fastq.gz ; do minimap2 -a -x sr reference.fasta $f 2> $(echo $f | cut -d'.' -f1).trimmed.pangolin_ref.minimap2.log | samtools view -bS - > $(echo $f | cut -d'.' -f1).trimmed.pangolin_ref.minimap2.bam ; done
```
