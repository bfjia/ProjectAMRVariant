## Table of contents
* [General info](#general-info)
* [Setup](#setup)
* [Usage](#usage)
* [Output](#output)
* [License](#license)

## General info
This repo contains the necessary workflow for detecting antimicrobial resistance gene (AMR) variant  in metagenomic data. 

Input: Paired-end illumina reads (.fastq, .fastq.gz), reference CARD database and SNP data (card.json, snp.txt).

Output: detected AMR SNPs


# Workflow:
1. Using card.json and snp.txt, a reference fasta file is created for each ARO that contain SNPs.
2. Reads are aligned to the reference fasta using bowtie2 in very-sensitive mode. 
3. Generate alignment pileup with samtools, include all positions. This step can be fine-tuned by specifying MAPQ quality(i.e. --MAPQ)
4. For each variant in snp.txt, scans the pileup file for that ARO+position, record the detected base, and depth. 
5. Calculate support values with base depth and count of variants detected.
6. Filter according to tuning parameters (i.e. --max-recall, --max-precision).
7. Output summary file(tab delimited).

## Setup

Dependencies: biopython=1.72, numpy, pandas, bcftools, samtools, and bowtie2.

It is suggested to run this program within a conda environment. These dependencies can be installed from `env.yaml` using `conda env create -f env.yaml`

Following environment installing, clone this repo into any directory and see `python main.py --help`

At somepoint, this will be build into a conda package. 

## Usage

`python main.py --forward /path/to/forward.fq --reverse /path/to/reverse.fq --card-json /path/to/card.json --card-snp path/to/snps.txt`

At somepoint, the default behavior for --card-json and --card-snp will be to download the latest copy from card.mcmaster.ca unless specified. For now, you must specify it's path.

# Example
See `test/` for examples. To run the workflow using these test files:

`python main.py --forward ./test/synthetic_1.fq --reverse ./test/synthetic_2.fq --temp ./example --card-json ./test/card.json --card-snp ./test/snps.txt --max-recall`

This will generate an output folder called `example` and a result file called `example_DetectedVariants.tsv`

## Output
The output is a tab delimitted file containing the following information:

| Column          | Description                        |
| --------------- | -----------------------------------|
| ARO     	      | ARO of resistance gene             |
| VariantClass    | CARD Model Type                    |
| VariantType     | Mutation Type                      |
| ResistantVariant| Flag for Resistant Variant         |
| SNP     	      | Detected SNP                       |
| Depth     	  | Total depth of the SNP position    |
| AbsSupport      | Number of reads matching the SNP   |
| RelativeSuppor  | %of resistant reads to total reads |
| INFO     	      | Raw counts of alternate variants   |

## Contribute

Any contributions are welcomed.
1. Fork it!
2. Create your feature branch: git checkout -b my-new-feature
3. Commit your changes: git commit -am 'Add some feature'
4. Push to the branch: git push origin my-new-feature
5. Submit a pull request :D


## License

GNU General Public License v3.0

Permissions of this strong copyleft license are conditioned on making available complete sourcecode of licensed works and modifications, which include larger works using a licensed work, under the same license. Copyright and license notices must be preserved. Contributors provide an express grant of patent rights.