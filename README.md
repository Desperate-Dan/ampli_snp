# primalscreen
A tool to find the minimal number of primer pairs required to cover a list of snps when given a primer bed file. Initially desgined to work for Monkeypox Virus amplicon based primer schemes, this will work with any amplicon based primer scheme suh as the ARTIC SARS-CoV-2 schemes.

This tool can take a json file of lineage defining mutations constructed like the Moneypox lineage.json (https://github.com/mpxv-lineages/lineage-designation/blob/master/auto-generated/lineages.json) and/or a csv of SNP positions and the lineage they define. It then compares this to a provided bed file of primer amplicon positions and returns a tsv of the primers that cover all SNPs in the CSV/json provided, as we as the minimal amplicon list to cover all SNPs for the lineages of interest, and finally the minimal list of amplicons that contian SNPs to differentiate all lineages of interest.

## Installation
```
git clone https://github.com/Desperate-Dan/primalscreen.git && cd primalscreen
pip install .
```

## Usage
```
usage: primalscreen [-h] [-c LINEAGE_CSV] [-j LINEAGE_JSON] [-l LIN_OF_INT] input_bed

Provide a an amplicon scheme bed file with either a json or csv of snps and you will recieve a tsv of amplicons that cover your SNPs of interest. If multiple
lineages are provided you can see the minimal amplicons that will allow distinguishing of the lineages based on these SNPs

positional arguments:
  input_bed             Input BED file to identify amplicons

optional arguments:
  -h, --help            show this help message and exit
  -c LINEAGE_CSV, --csv LINEAGE_CSV
                        CSV of SNPs of interest with format: SNP_position,lineage
  -j LINEAGE_JSON, --json LINEAGE_JSON
                        JSON of SNPs of interest in format of 'https://github.com/mpxv-lineages/lineage-designation/blob/master/auto-generated/lineages.json'
  -l LIN_OF_INT, --lineages LIN_OF_INT
                        Comma Separated list of lineages of interest from the input CSV or JSON. Default: all lineages in CSV and JSON                  
```
