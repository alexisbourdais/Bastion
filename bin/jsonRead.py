#!/usr/bin/env python

import json
import os

with open('short_summary.specific.pseudomonadales_odb10.TestBusco120.json', 'r') as file:
    with open("Test_lecture_json.csv", 'a') as results:
        data = json.load(file)

        #Lineage
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

        #Entetes
        print("Input_file\tDataset\tComplete\tSingle\tDuplicated\tFragmented\tMissing_n_markers\tScaffold N50\tContigs N50\tPercent gaps\tNumber of scaffolds", file=results)

        #Data
        print(f"{Input_file}\t{lineage}\t{Complete}\t{Single}\t{Duplicated}\t{Fragmented}\t{Missing_n_markers}\t{ScaffoldN50}\t{ContigsN50}\t{PercentGaps}\t{NumberScaffolds}", file=results)
