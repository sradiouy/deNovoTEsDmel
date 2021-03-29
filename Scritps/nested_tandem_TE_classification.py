#!/usr/bin/env python -W ignore::DeprecationWarning
# encoding: utf-8

## This script is to write down the nested and tandem-tes

import pandas as pd 
import subprocess
import sys
import os

def runcml(cml):
    """Ejecuta un comando"""
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    print(stderr)
    stdout = stdout.decode("utf-8")
    return stdout


def bedtoolsclosest(reference_gene_bed,te_bed,outfile,mode):
    """
    reference_gene_bed: es el archivo de genes de la referencia preparado para cada genoma
    te_bed son las coordenadas de los TE de cada genoma en el genoma de referencia
    """
    if mode == "overlap":
        commandLine = "bedtools closest -D ref -t all -a " + te_bed + " -b " + reference_gene_bed
    elif mode == "upstream": # -io ignoro overlap, ignoro -id downstream
        commandLine = "bedtools closest -D ref -t all -io -id -a " + te_bed + " -b " + reference_gene_bed
    else:
        # -io ignoro overlap, ignoro -iu upstream
        commandLine = "bedtools closest -D ref -t all -io -iu -a " + te_bed + " -b " + reference_gene_bed
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def parse_closest_files_te_te(closest_file):
    header = ["chr","start","end","TEID","len","strand","chr_1","start_1","end_1","TEID_1","len_1","strand_1","distance"]
    used_cols = ["chr","start","end","TEID","len","strand","TEID_1","distance"]
    df_bed = pd.read_table(closest_file,sep="\t",header=None,names=header)
    df_bed = df_bed[used_cols]
    return df_bed



def intersect_tes_tes_same_strain(te_annotation_bed,te_annotation_bed_2,outfile):
    """
    El objetivo de este paso es crear un archivo de interseccion entre los TEs de un mismo genoma, para poder determinar la existencia de nested TEs 
    """
    commandLine = "bedtools intersect -a " + te_annotation_bed + " -b " + te_annotation_bed_2 + " -wao "
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def parse_intersect_tes_tes_same_strain(intersection_file):
    header = ["chr","start","end","TEID","len","strand","chr_1","start_1","end_1","TEID_1","frame_1","strand_1","intersection"]
    used_cols = ["chr","start","end","TEID","len","strand","TEID_1","intersection"]
    df_bed = pd.read_table(intersection_file,sep="\t",header=None,names=header)
    df_bed = df_bed[used_cols]
    return df_bed


# Cargo el archivo log


logfile = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2]+ ".log"


### Primero que todo creo un directorio 

fo = open(logfile,"a")
log_line = "--- Determination of nested and tandem TE's ---\n"
fo.write(log_line)


te_original_strain_bed = sys.argv[3] 

fo = open(logfile,"a")
log_line = "--- Running bedtools intersect to obtain nested TE's ---\n"
fo.write(log_line)

te_intersect =intersect_tes_tes_same_strain(te_original_strain_bed,te_original_strain_bed,"te_intersect.txt")

fo = open(logfile,"a")
log_line = "--- Running bedtools intersect to obtain nested TE's: Done ---\n"
fo.write(log_line)


df_intersection_te_same_strain = parse_intersect_tes_tes_same_strain(te_intersect)

te_set = set(df_intersection_te_same_strain.TEID.tolist())


dict_te_len = {te:df_intersection_te_same_strain[df_intersection_te_same_strain["TEID"] == te].len.values[0] for te in te_set}

te_dict_full_nested = {te:[] for te in te_set}


for index,row in df_intersection_te_same_strain.iterrows():
    if row.TEID == row.TEID_1:
        nested = "NO"
    else:
        if row.len == row.intersection:
            te_dict_full_nested[row.TEID_1] += [row.TEID]

temp = {}
for key,values in te_dict_full_nested.items():
    if len(values) != 0:
        temp[key] = values

te_dict_full_nested = temp


partial_nested = []

for index,row in df_intersection_te_same_strain.iterrows():
    if row.TEID == row.TEID_1:
        pass
    else:
        if row.len != row.intersection:
            if row.TEID in te_dict_full_nested.keys():
                if not row.TEID_1 in te_dict_full_nested[row.TEID]:
                    if not row.TEID_1 in partial_nested:
                        partial_nested.append(row.TEID_1)
                    if not row.TEID in partial_nested:
                        partial_nested.append(row.TEID)
            else:
                if not row.TEID_1 in partial_nested:
                        partial_nested.append(row.TEID_1)
                if not row.TEID in partial_nested:
                    partial_nested.append(row.TEID)
                

nested_info = [] 
for key,values in te_dict_full_nested.items():
    for value in values:
        if value in partial_nested:
            nested_info.append([value,key,"Y"])
        else:
            nested_info.append([value,key,"N"])

nested_tes = [te[0] for te in nested_info]
for te in te_set:
    if te not in nested_tes:
        if te in partial_nested:
            nested_info.append([te,".","Y"])
        else:
            nested_info.append([te,".","N"])

df_nested = pd.DataFrame(nested_info,columns=["ID","inside_of","partial_nested"])

# Puede ser que este dentro de mas de un TE, en ese caso me quedo con uno
df_nested = df_nested.sort_values("ID")
df_nested = df_nested.drop_duplicates("ID",keep="first")


### Ahora veo informacion de tandem

fo = open(logfile,"a")
log_line = "--- Running bedtools closest to obtain the upstream closest TE for each TE (ignoring overlap and downstream) ---\n"
fo.write(log_line)

upstream = bedtoolsclosest(te_original_strain_bed,te_original_strain_bed,"upstream.txt","upstream")

fo = open(logfile,"a")
log_line = "--- Running bedtools closest to obtain the upstream closest TE for each TE (ignoring overlap and downstream): Done ---\n"
fo.write(log_line)

fo = open(logfile,"a")
log_line = "--- Running bedtools closest to obtain the downstream closest TE for each TE (ignoring overlap and upstream) ---\n"
fo.write(log_line)


downstream = bedtoolsclosest(te_original_strain_bed,te_original_strain_bed,"downstream.txt","downstream")

fo = open(logfile,"a")
log_line = "--- Running bedtools closest to obtain the downstream closest TE for each TE (ignoring overlap and upstream): Done ---\n"
fo.write(log_line)

df_upstream = parse_closest_files_te_te(upstream)
df_downstream = parse_closest_files_te_te(downstream)

## Analizo tandem en up

te_dict_up = {te:"" for te in te_set}

for index,row in df_upstream.iterrows():
    if row.TEID_1 == ".":
        pass
    else:
        distance = abs(row.distance)
        if distance <= 500:
            te_dict_up[row.TEID] = row.TEID_1 + ":" + str(distance)


### Elimino las key vacias

temp = {}
for key,values in te_dict_up.items():
    if len(values) != 0:
        temp[key] = values

te_dict_up = temp

### Hago lo mismo en down

te_dict_down = {te:"" for te in te_set}

for index,row in df_downstream.iterrows():
    if row.TEID_1 == ".":
        pass
    else:
        distance = abs(row.distance)
        if distance <= 500:
            te_dict_down[row.TEID] = row.TEID_1 + ":" + str(distance)


### Elimino las key vacias

temp = {}
for key,values in te_dict_down.items():
    if len(values) != 0:
        temp[key] = values

te_dict_down = temp


tandem_info = []
for te in te_set:
    tandem_up = "."
    tandem_down = "."
    if te in te_dict_down.keys():
        tandem_down =te_dict_down[te]
    if te in te_dict_up.keys():
        tandem_up =te_dict_up[te]
    tandem_info.append([te,tandem_up,tandem_down])

df_tandem = pd.DataFrame(tandem_info,columns=["ID","tandem_up","tandem_down"])


## Ahora junto la informacion de nested y tandem 

df_nested_tandem = df_nested.merge(df_tandem,on="ID")


nested_te_file = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_nested_tandem_info.tsv"


fo = open(logfile,"a")
log_line = "--- Writing nested and tandem information: " + nested_te_file + " ---\n"
fo.write(log_line)

df_nested_tandem.to_csv(nested_te_file,header=True,index=False,sep="\t")


os.remove("te_intersect.txt")
os.remove("upstream.txt")
os.remove("downstream.txt")

fo = open(logfile,"a")
log_line = "--- Writing nested and tandem information: Done ---\n"
fo.write(log_line)
