# CovidQuant
CovidQuant is a tool to quantify [PANGO lineage](https://cov-lineages.org/) probabilities from a mixed SARS-CoV-2 sequencing dataset (trimmed SAM/BAM) using the [PangoLEARN](https://github.com/cov-lineages/pangoLEARN) decision tree. CovidQuant loads the PangoLEARN decision tree, computes the number of reads that support each decision in the tree, and performs the following:

* Assign a sample to a single PANGO lineage (much like [Pangolin](https://github.com/cov-lineages/pangolin)), but directly using a trimmed SAM/BAM (rather than from a consensus genome sequence)
* Compute the overall probability of *every* PANGO lineage

## Installation
CovidQuant is written in Python 3 and depends on [pysam](https://pysam.readthedocs.io/). You can simply download [CovidQuant.py](CovidQuant.py) to your machine and make it executable:

```bash
wget "https://raw.githubusercontent.com/niemasd/CovidQuant/master/CovidQuant.py"
chmod a+x CovidQuant.py
sudo mv CovidQuant.py /usr/local/bin/CovidQuant.py # optional step to install globally
```

## Usage
CovidQuant can be used as follows:

```
usage: CovidQuant.py [-h] -i INPUT_ALIGNMENT -p PANGOLEARN_RULES [-o OUTPUT_ABUNDANCES] [--assign_top_lineage] [--output_decision_coverage OUTPUT_DECISION_COVERAGE] [-u]

optional arguments:
  -h, --help                                                     show this help message and exit
  -i INPUT_ALIGNMENT, --input_alignment INPUT_ALIGNMENT          Input Alignment File (SAM/BAM) (default: None)
  -p PANGOLEARN_RULES, --pangolearn_rules PANGOLEARN_RULES       Input PangoLEARN Decision Tree Rules (default: None)
  -o OUTPUT_ABUNDANCES, --output_abundances OUTPUT_ABUNDANCES    Output Abundances (TSV) (default: stdout)
  --assign_top_lineage                                           Optional: Assign top lineage using decision tree (will only be included in log output) (default: False)
  --output_decision_coverage OUTPUT_DECISION_COVERAGE            Optional: Output coverage of each PangoLEARN genome position (TSV) (default: None)
  -u, --update                                                   Update CovidQuant (default: False)
```

Note that the reads in the trimmed BAM/SAM must have been aligned to the exact reference genome used by Pangolin:

https://github.com/cov-lineages/pangolin/blob/master/pangolin/data/reference.fasta

If you have a trimmed BAM that you want to realign to the Pangolin reference genome, you can do so as follows:

1. Convert the BAM to FASTQ using [bedtools bamtofastq](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html):
  * `bedtools bamtofastq -i <INPUT_BAM> -fq <OUTPUT_FASTQ>`
2. Remap the FASTQ using [Minimap2](https://github.com/lh3/minimap2):
  * To output as a SAM (larger file): `minimap2 -a -x sr <REFERENCE_FASTA> <INPUT_FASTQ> > <OUTPUT_SAM>`
  * To output as a BAM (smaller file, needs [samtools](http://www.htslib.org/)): `minimap2 -a -x sr <REFERENCE_FASTA> <INPUT_FASTQ> | samtools view -bS - > <OUTPUT_BAM>`

# Citing CovidQuant
If you use CovidQuant in your work, please cite this GitHub repo (paper TBD).
