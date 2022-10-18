#!/usr/bin/env python3

# import some bits
import sys
import os
import argparse
from collections import defaultdict

import function_file as funk

def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(description='Provide a an amplicon scheme bed file with either a json or csv of snps and you will recieve a tsv of amplicons that cover your SNPs of interest. If multiple lineages are provided you can see the minimal amplicons that will allow distinguishing of the lineages based on these SNPs')

    parser.add_argument('input_bed', help = "Input BED file to identify amplicons", action = 'store')
    parser.add_argument('-c', '--csv', help = "CSV of SNPs of interest with format: SNP_position,lineage", action = 'store', dest = 'lineage_csv')
    parser.add_argument('-j', '--json', help = "JSON of SNPs of interest in format of 'https://github.com/mpxv-lineages/lineage-designation/blob/master/auto-generated/lineages.json'", action = 'store', dest = 'lineage_json')
    parser.add_argument('-l', '--lineages', help = "Comma Separated list of lineages of interest from the input CSV or JSON. Default: all lineages in CSV and JSON", action = 'store', dest = 'lin_of_int', type = str, default = set())

    if len(sysargs) < 1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)


    #Here is the business

    output_amplicon_dict = funk.bed_file_reader(args.input_bed)

    output_def_snps_csv, output_snp_lin_csv = funk.lineage_csv_parser(args.lineage_csv)

    output_def_snps_json, output_snp_lin_json = funk.lineage_json_parser(args.lineage_json)

    funk.amplicons_vs_snps(output_amplicon_dict, output_snp_lin_csv, output_snp_lin_json, args.lin_of_int, output_def_snps_json, output_def_snps_csv)
    

if __name__ == '__main__':
    main()