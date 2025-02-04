#!/usr/bin/env python

import argparse
import csv
import os

#######################
### Argument parser ###
#######################
ap = argparse.ArgumentParser()

### Files
ap.add_argument("--busco", default='', required=False)
ap.add_argument("--quast", default='', required=False)
ap.add_argument("--eukcc", default='', required=False)
ap.add_argument("--gunc", default='', required=False)
ap.add_argument("--checkm1", default='', required=False)
ap.add_argument("--checkm2", default='', required=False)
ap.add_argument("--omark", default='', required=False)
ap.add_argument("--kmerfinder", default='', required=False)
ap.add_argument("--gtdbtk", default='', required=False)
ap.add_argument("--physeter", default='', required=False)
ap.add_argument("--kraken", default='', required=False)
ap.add_argument("--plasmidfinder", default='', required=False)

args = vars(ap.parse_args())

#################
### Variables ###
#################
genome_list = []

checkm1_complet_dic = {}
checkm1_contam_dic = {}
checkm1_str_dic = {}

gunc_css_dic = {}
gunc_contam_dic = {}
gunc_RRS_dic = {}
gunc_status_dic = {}

busco_placement_dic = {}
busco_complet_dic = {}
busco_single_dic = {}
busco_duplicated_dic = {}
busco_fragmented_dic = {}
busco_missing_dic = {}

quast_contig_dic = {}
quast_length_dic = {}
quast_gc_dic = {}
quast_n50_dic = {}

checkm2_complet_dic = {}
checkm2_contam_dic = {}

gtdbtk_placement_dic = {}
gtdbtk_reference_dic = {}

eukcc_complet_dic = {}
eukcc_contam_dic = {}

omark_placement_dic = {}
omark_score_dic = {}

physeter_contam_dic = {}
physeter_placement_dic = {}

kraken_contam_dic = {}
kraken_placement_dic = {}

kmerfinder_RefSeq_id_dic = {}
kmerfinder_taxonomy_dic = {}
kmerfinder_tot_template_Coverage_dic = {}
kmerfinder_tot_query_Coverage_dic = {}

omark_ref_dic = {}
omark_score_dic = {}
omark_contam_dic = {}

plasmidfinder_hit_dic = {}

##############
### Script ###
##############

#Checkm1
if os.path.isfile(args['checkm1']):
    checkm1 = csv.reader(open(args['checkm1'], "r"), delimiter=',')
    for line in checkm1:
        genome_id = line[0]
        completeness = line[1]
        contam = line[2]
        Str = line[3]

        genome_list.append(genome_id)
        checkm1_complet_dic[genome_id] = completeness
        checkm1_contam_dic[genome_id] = contam
        checkm1_str_dic[genome_id] = Str

#Gunc
if os.path.isfile(args['gunc']):
    gunc = csv.reader(open(args['gunc'], "r"), delimiter='\t')
    for line in gunc:
        if ('n_genes_called' in line):
            continue
        else:
            genome_id = line[0].replace(".", "")
            css = line[7]
            contam = line[8]
            rrs = line[11]
            status = line[12]

            gunc_contam_dic[genome_id] = contam
            gunc_RRS_dic[genome_id] = rrs
            gunc_css_dic[genome_id] = css 
            gunc_status_dic[genome_id] = status 

#Busco
if os.path.isfile(args['busco']):
    busco = csv.reader(open(args['busco'], "r"), delimiter='\t')
    for line in busco:
        if ('Input_file' in line):
            continue
        else:
            genome_id = line[0].split('.')[0]
            placement = line[1]
            completness = line[2]
            single = line[3]
            duplicated = line[4]
            fragmented = line[5]
            missing = line[6]

            busco_placement_dic[genome_id] = placement
            busco_complet_dic[genome_id] = completness
            busco_single_dic[genome_id] = single
            busco_duplicated_dic[genome_id] = duplicated
            busco_fragmented_dic[genome_id] = fragmented
            busco_missing_dic[genome_id] = missing

#Quast
if os.path.isfile(args['quast']):
    quast = csv.reader(open(args['quast'], "r"), delimiter='\t')
    for line in quast:
        if ('Assembly' in line):
            continue
        else:
            genome_id = line[0]
            contigs = line[1]
            length = line[7]
            gc = line[16]
            n50 = line[17]

            quast_contig_dic[genome_id] = contigs
            quast_length_dic[genome_id] = length
            quast_gc_dic[genome_id] = gc
            quast_n50_dic[genome_id] = n50

#Checkm2
if os.path.isfile(args['checkm2']):
    checkm2 = csv.reader(open(args['checkm2'], "r"), delimiter='\t')
    for line in checkm2:
        if ('Name' in line):
            continue
        else:
            genome_id = line[0]
            completeness = line[1]
            contam = line[2]

            checkm2_complet_dic[genome_id] = completeness
            checkm2_contam_dic[genome_id] = contam

#Eukcc
if os.path.isfile(args['eukcc']):
    eukcc = csv.reader(open(args['eukcc'], "r"), delimiter='\t')
    for line in eukcc:
        if ('bin' in line):
            continue
        else:
            genome_id = line[0].split('.')[0]
            complet = line[1]
            contam = line[2]

            eukcc_complet_dic[genome_id] = complet
            eukcc_contam_dic[genome_id] = contam

#GTDBTK
if os.path.isfile(args['gtdbtk']):
    gtdbtk = csv.reader(open(args['gtdbtk'], "r"), delimiter='\t')
    for line in gtdbtk:
        if ('user_genome' in line):
            continue
        else:
            genome_id = line[0]
            placement = line[1]
            reference = line[2]

            gtdbtk_placement_dic[genome_id] = placement
            gtdbtk_reference_dic[genome_id] = reference

#Physeter
if os.path.isfile(args['physeter']):
    physeter = csv.reader(open(args['physeter'], "r"), delimiter='\t')
    for line in physeter:
        genome_id = line[0].replace("-abbr-split-kraken", "")
        placement = line[1]
        contam = line[3]
        physeter_contam_dic[genome_id] = contam
        physeter_placement_dic[genome_id] = placement

#Kraken
if os.path.isfile(args['kraken']):
    kraken = csv.reader(open(args['kraken'], "r"), delimiter='\t')
    for line in kraken:
        genome_id = line[0].replace("-split.report", "")
        placement = line[1]
        contam = line[3]
        kraken_contam_dic[genome_id] = contam
        kraken_placement_dic[genome_id] = placement   

#Kmerfinder
if os.path.isfile(args['kmerfinder']):
    kmerfinder = csv.reader(open(args['kmerfinder'], "r"), delimiter='\t')
    for line in kmerfinder:
        if ('# Assembly' in line):
            continue
        else:
            genome_id = line[0]
            RefSeq_id = line[1]
            taxonomy = line[17]
            tot_template_Coverage = line[10]
            tot_query_coverage = line[9]
        
            kmerfinder_RefSeq_id_dic[genome_id] = RefSeq_id
            kmerfinder_taxonomy_dic[genome_id] = taxonomy
            kmerfinder_tot_template_Coverage_dic[genome_id] = tot_template_Coverage
            kmerfinder_tot_query_Coverage_dic = tot_query_coverage

#Omark
if os.path.isfile(args['omark']):
    omark = csv.reader(open(args['omark'], "r"), delimiter='\t')
    for line in omark:
        genome_id = line[0]
        placement = line[1]
        score = line[2]
        contam = line[3]

        omark_ref_dic[genome_id] = placement
        omark_score_dic[genome_id] = score
        omark_contam_dic[genome_id] = contam

#Plasmidfinder
if os.path.isfile(args['plasmidfinder']):
    plasmidfiner = csv.reader(open(args['plasmidfinder'], "r"), delimiter='\t')
    for line in plasmidfiner:
        genome_id = line[0]
        plasmid_hit = line[1]

        plasmidfinder_hit_dic[genome_id] = plasmid_hit

### Results
results_file=open("Bastion_FinalReport.tsv", "a")
print("Genome\t\
Checkm_completeness\tCheckm_contamination\tCheckm_str\t\
Busco_placemenent\tBusco_completeness\tBusco_duplicate\t\
GUNC_CSS\tGUNC_RRS\tGUNC_status\tGUNC_conta\t\
Physeter_placement\tPhyseter_contamination\t\
Kraken_placement\tKraken_contamination\t\
Checkm2_completeness\tCheckm2_contamination\t\
Gtdbtk_placement\tGtdbtk_reference\t\
Eukcc2_contamination\tEukcc2_completeness\t\
quast_#contigs\tquast_tot_length\tquast_GC\tquast_N50\t\
kmerfinder_RefSeq_id\tkmerfinder_taxonomy\tKmerfinder_Tot_Template_Coverage\ttkmerfinder_Tot_Query_Coverage\t\
Omark_Main_species\tOmark_score\tOmark_Contam\t\
Plasmidfinder", file=results_file)

for genome in genome_list:
    checkm1_complet=checkm1_complet_dic.get(genome)
    checkm1_contam=checkm1_contam_dic.get(genome)
    checkm1_str=checkm1_str_dic.get(genome)

    gunc_css=gunc_css_dic.get(genome)
    gunc_contam=gunc_contam_dic.get(genome)
    gunc_RRS=gunc_RRS_dic.get(genome)
    gunc_status=gunc_status_dic.get(genome)
    
    busco_placement=busco_placement_dic.get(genome)
    busco_complet=busco_complet_dic.get(genome)
    busco_single=busco_single_dic.get(genome)
    busco_duplicated=busco_duplicated_dic.get(genome)
    busco_fragmented=busco_fragmented_dic.get(genome)
    busco_missing=busco_missing_dic.get(genome)

    quast_contig=quast_contig_dic.get(genome)
    quast_length=quast_length_dic.get(genome)
    quast_gc=quast_gc_dic.get(genome)
    quast_n50=quast_n50_dic.get(genome)

    checkm2_complet=checkm2_complet_dic.get(genome)
    checkm2_contam=checkm2_contam_dic.get(genome)

    gtdbtk_placement=gtdbtk_placement_dic.get(genome)
    gtdbtk_reference=gtdbtk_reference_dic.get(genome)

    eukcc_complet=eukcc_complet_dic.get(genome)
    eukcc_contam=eukcc_contam_dic.get(genome)

    physeter_contam=physeter_contam_dic.get(genome)
    physeter_placement=physeter_placement_dic.get(genome)

    kraken_contam=kraken_contam_dic.get(genome)
    kraken_placement=kraken_placement_dic.get(genome)

    kmerfinder_RefSeq_id=kmerfinder_RefSeq_id_dic.get(genome)
    kmerfinder_taxonomy=kmerfinder_taxonomy_dic.get(genome)
    kmerfinder_tot_template_Coverage=kmerfinder_tot_template_Coverage_dic.get(genome)
    kmerfinder_tot_query_Coverage=kmerfinder_tot_query_Coverage_dic.get(genome)

    omark_ref=omark_ref_dic.get(genome)
    omark_score=omark_score_dic.get(genome)
    omark_contam=omark_contam_dic.get(genome)

    plasmid_results=plasmidfinder_hit_dic.get(genome)

    print(f"{genome}\t\
{checkm1_complet}\t{checkm1_contam}\t{checkm1_str}\t\
{busco_placement}\t{busco_complet}\t{busco_duplicated}\t\
{gunc_css}\t{gunc_RRS}\t{gunc_status}\t{gunc_contam}\t\
{physeter_placement}\t{physeter_contam}\t\
{kraken_placement}\t{kraken_contam}\t\
{checkm2_complet}\t{checkm2_contam}\t\
{gtdbtk_placement}\t{gtdbtk_reference}\t\
{eukcc_contam}\t{eukcc_complet}\t\
{quast_contig}\t{quast_length}\t{quast_gc}\t{quast_n50}\t\
{kmerfinder_RefSeq_id}\t{kmerfinder_taxonomy}\t{kmerfinder_tot_template_Coverage}\t\{kmerfinder_tot_query_Coverage}t\
{omark_ref}\t{omark_score}\t{omark_contam}\t\
{plasmid_results}", file=results_file)