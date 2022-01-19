# AMR Variant Detection In Metagenomic Data

Variant detection workflow for AMR genes in metagenomic data

Entry point: main.py
Input: Paired-end illumina reads (.fastq, .fastq.gz), reference CARD database and SNP data (card.json, snp.txt).
Output: detected AMR SNPs


Workflow:
1. Using card.json and snp.txt, a reference fasta file is created for each ARO that contain SNPs.
2. Reads are aligned to the reference fasta using bowtie2 in very-sensitive mode. 
3. Generate alignment pileup with samtools, include all positions. This step can be fine-tuned by specifying MAPQ quality(i.e. --MAPQ)
4. For each variant in snp.txt, scans the pileup file for that ARO+position, record the detected base, and depth. 
5. Calculate support values with base depth and count of variants detected.
6. Filter according to tuning parameters (i.e. --max-recall, --max-precision).
7. Output summary file(tab delimited).

## Installation

Dependencies: biopython=1.72, numpy, pandas, bcftools, samtools, and bowtie2.

It is suggested to run this program within a conda environment. These dependencies can be installed from `env.yaml` using `conda env create -f env.yaml`

Following environment installing, clone this repo into any directory and see `python main.py --help`

At somepoint, this will be build into a conda package. 

## Example Execution

`python main.py --forward /path/to/forward.fq --reverse /path/to/reverse.fq --card-json /path/to/card.json --card-snp path/to/snps.txt`

At somepoint, the default behavior for --card-json and --card-snp will be to download the latest copy from card.mcmaster.ca unless specified. For now, you must specify it's path.

See `test/` for examples. To run the workflow using these test files:

`python main.py --forward ./test/synthetic_1.fq --reverse ./test/synthetic_2.fq --temp ./example --card-json ./test/card.json --card-snp ./test/snps.txt --max-recall`

This will generate an output folder called `example` and a result file called `example_DetectedVariants.tsv`

Columns of output tsv
ARO	- ARO of resistance gene
VariantClass - CARD Model Type
VariantType	- Mutation Type
ResistantVariant - Flag for Resistant Variant
SNP	- Detected SNP
Depth - Total depth of the SNP position
AbsSupport - Number of reads matching the SNP
RelativeSupport	- %of resistant reads to total reads
INFO - Raw counts of alternate variants. 