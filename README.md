# CovidQuant
CovidQuant is a tool to quantify [PANGO lineage](https://cov-lineages.org/) probabilities from a mixed SARS-CoV-2 sequencing dataset (trimmed SAM/BAM) using the [PangoLEARN](https://github.com/cov-lineages/pangoLEARN) decision tree. CovidQuant loads the PangoLEARN decision tree, computes the number of reads that support each decision in the tree, and performs the following:

* Compute the number of reads that support an `A`, `C`, `G`, `T`, or `-` at every reference genome position utilized by PangoLEARN
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

## Approach
The PangoLEARN decision tree has internal nodes labeled by `(genome position, symbol)` pairs (where `symbol` is `A`, `C`, `G`, `T`, or `-`), and it has leaves labeled by PANGO lineage (each leaf is labeled by a single PANGO lineage, but a single PANGO lineage can label multiple leaves). The general idea behind CovidQuant is to estimate support for each decision in the PangoLEARN decision tree directly using reads mapped to the reference genome (rather than from a consensus sequence). For each decision, CovidQuant counts the number of reads that support vs. contradict that decision and estimates an empirical probability.

### Single PANGO Lineage Assignment
Starting at the root, CovidQuant iteratively takes the highest-probability decision until it reaches a leaf, and it outputs the PANGO lineage labeling that leaf as the assigned lineage.

### Probability of Every PANGO Lineage
Each individual decision (which is represented as a single edge in the decision tree) has an empirical probability: the number of reads supporting that decision divided by the number of reads that either support or contradict that decision. The probability of a path from the root to a leaf is the product of the probabilities of the edges along that path. CovidQuant defines the probability of a single PANGO lineage *x* as the sum of the probabilities of all root-to-leaf paths in which the leaf is labeled by *x*.

# Citing CovidQuant
If you use CovidQuant in your work, please cite this GitHub repo (paper TBD).
