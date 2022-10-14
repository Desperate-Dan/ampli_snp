from Bio import SeqIO
import re
from collections import defaultdict
import json

test_bed = "MPXV.primer.bed"
lineage_def = "lineages.json"
lineage_csv = "MPX_reconstruction_with_only_def.csv"
 

def bed_file_reader(input_bed):
    #Read a bed file in format CHR, start_pos, end_pos, primer_id, etc and return an amplicon dictionary
    #Currently does not handle "alt" primers
    bed = open(input_bed)
    primer_list = []
    amplicon_dict = defaultdict(list)
    for line in bed:
        #Take each line of the bed file, strip \n, split on tab, keep first 4 elements to deal with 4 or 6 element bed files
        primer_list.append(line.rstrip("\n").split("\t")[:4])
    #Create a dictionary of amplicons, taking the 3' coordinates of each primer pair.    
    for primer in primer_list:
        if re.match(".*_LEFT", primer[-1]):
            amplicon_dict[primer[-1].rstrip("_LEFT")].append(int(primer[2]))
        elif re.match(".*_RIGHT", primer[-1]):
            amplicon_dict[primer[-1].rstrip("_RIGHT")].append(int(primer[1]))
    
    return (amplicon_dict)
 

def lineage_json_parser(lineage_json):
    #Parse the lineage json from "https://github.com/mpxv-lineages/lineage-designation/blob/master/auto-generated/lineages.json"
    data = json.load(open(lineage_json))
    defining_snps_dict_json = defaultdict(list)
    #Loop through the different lineages in the json...
    for lineage in data["lineages"]:
        #...and add their defining positions to a new dictionary.
        for def_snp in lineage["defining_snps"]:
            defining_snps_dict_json[lineage["name"]].append(def_snp["pos"])

    return defining_snps_dict_json

def lineage_csv_parser(lineage_csv):
    #Parse an input csv of mutations - be careful though, positions are hard coded...
    data = open(lineage_csv)
    defining_snps_dict_csv = defaultdict(list)
    #Skip the header
    next(data, None)
    for line in data:
        line = line.rstrip("\n").split(",")
        defining_snps_dict_csv[line[3]].append(int(line[2]))
    
    return defining_snps_dict_csv
    
def snp_lineage_inverter(lineage_json, lineage_csv):
    #Want to create a new dictionary with SNP keys and lineage values.
    #This could definitely be incorporated into the two snp defining functions above.
    data = open(lineage_csv)
    snp_lineage_dict = {}
    next(data, None)
    for line in data:
        line = line.rstrip("\n").split(",")
        snp_lineage_dict[int(line[2])] = line[3]
    data = json.load(open(lineage_json))
    defining_snps_dict_json = defaultdict(list)
    #Loop through the different lineages in the json...
    for lineage in data["lineages"]:
        #...and add their defining positions to a new dictionary.
        for def_snp in lineage["defining_snps"]:
            snp_lineage_dict[def_snp["pos"]] = lineage["name"]
            #defining_snps_dict_json[lineage["name"]].append(def_snp["pos"])
    return snp_lineage_dict

def amplicons_vs_snps(amplicon_dict, snp_lineage_dict,
                      defining_snps_dict_json = defaultdict(list), defining_snps_dict_csv = defaultdict(list)):
    #Combining the dictionaries made from the json and the csv - might be defferential to y values ( x | y ) so need to be careful
    defining_snps_dict = (defining_snps_dict_json | defining_snps_dict_csv)
    #print(defining_snps_dict)
    #Take the dict of snp positions and loop through the lineages.
    lineage_primer_dict = defaultdict(dict)
    all_primer_list = []
    for lineage in defining_snps_dict:
        #For each SNP set up a new dict to contain primers.
        snp_primer_dict = defaultdict(list)
        for snp in defining_snps_dict[lineage]:
            #If the SNP falls within an amplicon, add the primer pair to the lineage_primer_dict.
            lineage_primer_dict[lineage] = snp_primer_dict
            for primer in amplicon_dict:
                #If the "amplicon" does not have two coordinates, skip it (checks for bed file inconsistencies)
                if len(amplicon_dict[primer]) == 2:
                    if (snp > amplicon_dict[primer][0]) & (snp < amplicon_dict[primer][1]):
                        snp_primer_dict[snp].append(primer)
                        #extra list containing all primers to use for counts.
                        all_primer_list.append(primer)
                else:
                    continue

    #Write tsv of all amplicons that cover the input SNPs
    output_tsv = open("primer_pair_list_NG.tsv", "w")
    output_tsv.write("Lineage\tSNP_pos\tPrimer_Pairs\n")    
    for lineage in lineage_primer_dict:
        for snp in lineage_primer_dict[lineage]:
            output_tsv.write("%s\t%s\t%s\n" % (lineage,snp,lineage_primer_dict[lineage][snp]))
    output_tsv.close()
    
    #Create the minimal amplicon list
    #The idea here is to create a list of all primers needed, get counts for each primer for each SNP compared to it, keep the most common, then unique the list.
    #set up container dictionaries and lists
    primer_counts = {}
    snp_primer_counts = defaultdict(list)
    snp_primer_index = {}
    full_snp_primer_dict = defaultdict(list)
    minimal_primer_list = []
    amplicon_lineage_dict = defaultdict(set)
    #This should eventually be user defined but currently is all lineages with defining snps
    lineages_of_interest = set()
    
    #Make the overly complicated dictionaries and lists...
    for lineage in lineage_primer_dict:
        lineages_of_interest.add(lineage)
        for snp in lineage_primer_dict[lineage]:
            for i,primer in enumerate(lineage_primer_dict[lineage][snp]):
                #Count occurances of each primer in full list.
                primer_counts[primer] = all_primer_list.count(primer)
                #Then create a dict with a count for each primer with consistent index position as snp - primer list
                snp_primer_counts[snp].append(primer_counts[lineage_primer_dict[lineage][snp][i]])
                #Create the dict with consistent primer indexing
                full_snp_primer_dict[snp].append(primer)
                #Get the lineages that a particular amplicon covers
                amplicon_lineage_dict[primer].add(lineage)
            #Get the index of the primer that appears most often per snp    
            snp_primer_index[snp] = snp_primer_counts[snp].index(max(snp_primer_counts[snp]))
            #Append a list of the most often seen primer per snp
            minimal_primer_list.append(full_snp_primer_dict[snp][snp_primer_index[snp]])
    #"Unique" this list to get the minimal primer set
    #print("This is the minimal list of primers to cover every snp: %s" % set(minimal_primer_list))        
    
    #print(snp_primer_counts)
    #print(full_snp_primer_dict)
   
    #Should there be a check for potential "pool A and B" conflicts?
                
    #print (amplicon_lineage_dict)
    #print (lineages_of_interest)
    
    #Modified this from https://stackoverflow.com/questions/9281788/get-longest-element-in-dict/9281968#9281968
    amplicon_lineage_dict = (sorted(((len(v), v, k) for k, v in amplicon_lineage_dict.items()),reverse = True))
    
    
    least_primer_lineage_list = []
    covered_lineage_set = set()
    for i in amplicon_lineage_dict:
        if (covered_lineage_set == lineages_of_interest):
            break
        elif (covered_lineage_set.issuperset(i[1])):
            continue
        else:
            covered_lineage_set = covered_lineage_set.union(i[1])
            least_primer_lineage_list.append(i[2])
    
    #print (least_primer_lineage_list)
    
    #print (full_snp_primer_dict)
    #print("\n")
    #print (snp_primer_index)
    #print("\n")
    #print (lineage_primer_dict)
    #print("\n")
    #print (minimal_primer_list)
    
    return output_tsv, defining_snps_dict, least_primer_lineage_list
        
output_amplicon_dict = bed_file_reader(test_bed)

output_def_snps_csv = lineage_csv_parser(lineage_csv)

output_def_snps_json = lineage_json_parser(lineage_def)

output_snp_lin = snp_lineage_inverter(lineage_def, lineage_csv)

a, b, c = amplicons_vs_snps(output_amplicon_dict, output_snp_lin, output_def_snps_json, output_def_snps_csv)
