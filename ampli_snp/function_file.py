from Bio import SeqIO
import re
from collections import defaultdict
import json

#test_bed = "MPXV.primer.bed"
#lineage_def = "lineages.json"
#lineage_csv = "MPX_reconstruction_def.csv"



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
    defining_snps_dict_json = defaultdict(list)
    snp_lineage_dict_json = {}
    if lineage_json:
        data = json.load(open(lineage_json))
        #Loop through the different lineages in the json...
        for lineage in data["lineages"]:
            #...and add their defining positions to a new dictionary.
            for def_snp in lineage["defining_snps"]:
                defining_snps_dict_json[lineage["name"]].append(def_snp["pos"])
                snp_lineage_dict_json[def_snp["pos"]] = lineage["name"]

    return defining_snps_dict_json, snp_lineage_dict_json



def lineage_csv_parser(lineage_csv):
    #Parse an input csv of mutations - be careful though, positions are hard coded...
    defining_snps_dict_csv = defaultdict(list)
    snp_lineage_dict_csv = {}
    if lineage_csv:
        data = open(lineage_csv)    
        #Skip the header
        next(data, None)
        for line in data:
            line = line.rstrip("\n").split(",")
            defining_snps_dict_csv[line[1]].append(int(line[0]))
            snp_lineage_dict_csv[int(line[0])] = line[1]
    
    return defining_snps_dict_csv, snp_lineage_dict_csv



def amplicons_vs_snps(amplicon_dict, snp_lineage_dict_csv, snp_lineage_dict_json, lineages_of_interest,
                      defining_snps_dict_json = defaultdict(list), defining_snps_dict_csv = defaultdict(list)):
    #Combining the dictionaries made from the json and the csv - might be defferential to y values ( x | y ) so need to be careful
    defining_snps_dict = {**defining_snps_dict_json, **defining_snps_dict_csv}
    snp_lineage_dict = {**snp_lineage_dict_csv, **snp_lineage_dict_json}
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
    output_tsv = open("primer_pair_list.tsv", "w")
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
    
    #Parse the input lineages or select all lineages with defining snps if none given
    if lineages_of_interest:
        print(lineages_of_interest)
        lineages_of_interest = set(lineages_of_interest.split(","))
        print(lineages_of_interest)
        #for lineage in lineages_of_interest:

    else:
        for lineage in lineage_primer_dict:
            lineages_of_interest.add(lineage)
    
    #Make the overly complicated dictionaries and lists based on the lineages of interest
    for lineage in lineages_of_interest:
        for snp in lineage_primer_dict[lineage]:
            for i,primer in enumerate(lineage_primer_dict[lineage][snp]):
                #Count occurances of each primer in full list
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
    print("This is the minimal list of primers to cover every snp: %s" % set(minimal_primer_list))        

    
    #Next section is to generate the minimal primer set that captures at least one snp from all given lineages
    #Broadly this sorts the dictionary of sets of lineages covered by an amplicon from most to least
    #It then adds non-redundant primers to a list until all lineages are represented, checking for any "pool conflicts" 
    
    #Modified this from https://stackoverflow.com/questions/9281788/get-longest-element-in-dict/9281968#9281968
    #This now sorts the diictionary by the length of the values (note it's no longer actually a dictionary)
    amplicon_lineage_dict = (sorted(((len(v), v, k) for k, v in amplicon_lineage_dict.items()),reverse = True))
    
    #Create the lists and sets required.
    least_primer_lineage_list = []
    covered_lineage_set = set()
    for i in amplicon_lineage_dict:
        #First check if for any primers that cover all the lineages straight away
        if lineages_of_interest.issubset(i[1]):
            least_primer_lineage_list.append(i[2])
            break
        else:
            #If no easy match move on... Check if the covered lineages match with the ones we want to cover, when they do job's a good'un        
            if (covered_lineage_set.issuperset(lineages_of_interest)):
                break
            #If the currently covered lineages are a superset of the newly considered amplicon set, ignore it and continue
            elif (covered_lineage_set.issuperset(i[1])):
                continue
            #If it's not redundant to add it check if there are any primer pool conflicts
            else:
                #Check if the new amplicon is compatible with those already in the set
                if covered_lineage_set:
                    amp_no = int(i[2].split("_")[1])
                    if (re.match((".*" + str(amp_no + 1) + "'.*"), str(least_primer_lineage_list))):
                        print(least_primer_lineage_list)
                        print("\nAmplicon pool conflict found for %s and %s" % (i[2],amp_no + 1))
                        continue
                    elif (re.match((".*" + str(amp_no - 1) + "'.*"), str(least_primer_lineage_list))):
                        print("Amplicon pool conflict found for %s and %s" % (i[2],amp_no - 1))
                        continue                
                    else:
                        covered_lineage_set = covered_lineage_set.union(i[1])
                        least_primer_lineage_list.append(i[2])
                #This adds the first amplicon to allow subsequent checks
                else:    
                    covered_lineage_set = covered_lineage_set.union(i[1])
                    least_primer_lineage_list.append(i[2])
    
    print ("\nThis is the minimal list of primers to cover every lineage of interest: %s\n" % least_primer_lineage_list)
    
    
    return output_tsv, defining_snps_dict, least_primer_lineage_list