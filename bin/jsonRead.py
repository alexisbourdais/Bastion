#!/usr/bin/env python

###############
### Modules ###
###############
import json
import os
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

Input_file=(os.path.basename(data['parameters']['in']))
lineage=(data['lineage_dataset']['name'])
Complete=(data['results']['Complete percentage'])
Single=(data['results']['Single copy percentage'])
Duplicated=(data['results']['Multi copy percentage'])
Fragmented=(data['results']['Fragmented percentage'])
Missing_n_markers=(data['results']['Missing percentage']) 
ScaffoldN50=(data['metrics']['Scaffold N50'])
ContigsN50=(data['metrics']['Contigs N50'])
PercentGaps=(data['metrics']['Percent gaps']).replace("%", "")
NumberScaffolds=(data['metrics']['Number of scaffolds'])
Scores_archaea_odb10=""
Scores_bacteria_odb10=""
Scores_eukaryota_odb10=""

results_file=open(f"{Input_file}".replace(f".{args['format']}", "_BuscoResults.csv"), "a")

#Entetes
print("Input_file\tDataset\tComplete\tSingle\tDuplicated\tFragmented\tMissing_n_markers\tScaffold N50\tContigs N50\tPercent gaps\tNumber of scaffolds", file=results_file)

#Data
print(f"{Input_file}\t{lineage}\t{Complete}\t{Single}\t{Duplicated}\t{Fragmented}\t{Missing_n_markers}\t{ScaffoldN50}\t{ContigsN50}\t{PercentGaps}\t{NumberScaffolds}", file=results_file)