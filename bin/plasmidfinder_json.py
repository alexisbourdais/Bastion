#!/usr/bin/env python

###############
### Modules ###
###############

import json
import argparse

#######################
### Argument parser ###
#######################
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--input", required=True, help="Input file")
ap.add_argument("-f", "--format", required=True, help="Format file")

args = vars(ap.parse_args())

##############
### Script ###
##############
data = json.load(open(args['input'], 'r'))

Input_file=(data['plasmidfinder']['user_input']['filename(s)'][0])
name=f"{Input_file}".replace(f".{args['format']}", "")

results_entero=(data['plasmidfinder']['results']['Enterobacteriales']['enterobacteriales'])
results_Inc18=(data['plasmidfinder']['results']['Gram Positive']['Inc18'])
results_NT_Rep=(data['plasmidfinder']['results']['Gram Positive']['NT_Rep'])
results_NT_Rep1=(data['plasmidfinder']['results']['Gram Positive']['Rep1'])
results_NT_Rep2=(data['plasmidfinder']['results']['Gram Positive']['Rep2'])
results_NT_Rep3=(data['plasmidfinder']['results']['Gram Positive']['Rep3'])
results_NT_RepA_N=(data['plasmidfinder']['results']['Gram Positive']['RepA_N'])
results_NT_RepL=(data['plasmidfinder']['results']['Gram Positive']['RepL'])
results_NT_Rep_trans=(data['plasmidfinder']['results']['Gram Positive']['Rep_trans'])

if results_entero == "No hit found" and results_Inc18 == "No hit found" and results_NT_Rep == "No hit found" and results_NT_Rep1 == "No hit found" and results_NT_Rep2 == "No hit found" and results_NT_Rep3 == "No hit found" and results_NT_RepA_N == "No hit found" and results_NT_RepL == "No hit found" and results_NT_Rep_trans == "No hit found":
    plasmid="No hit"
else:
    plasmid="Potential plasmid"

results_file=open(f"{name}_plasmidfinder_results.tsv", "a")

#Data
print(f"{name}\t{plasmid}", file=results_file)