#!/usr/bin/env python -W ignore::DeprecationWarning
# encoding: utf-8

# Una vez obtenido los fasta de TEs a utilizar necesitamos encontrarlos en los genomas de referencia. Para eso utilizaremos minimap2 en la modalidad splice.

##############################################################################################################################
##############################################################################################################################

import subprocess
import pandas as pd
import pysam
from Bio import SeqIO
import numpy as np 
import sys
import os 
from time import gmtime, strftime
from collections import defaultdict

def ObtainSeq(genome):
    """
    Transforme a genome in fasta format into a dictionary..
    :param: file: genome: genome file in fasta format
    :return: genome sequence dictionary
    """
    handle = open(genome, "rU")
    record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()
    return record_dict

def get_m_values(cigartuples):
    m_values = []
    suma = 0
    for index,cigar_tuple in enumerate(cigartuples):
        if cigar_tuple[0] == 0 or cigar_tuple[0] == 1 or cigar_tuple[0] == 4:
            suma += cigar_tuple[1]
        elif cigar_tuple[0] == 2:
            pass
        elif cigar_tuple[0] == 3:
            m_values.append(suma)
            suma = 0
        if index == len(cigartuples)-1:
            m_values.append(suma)
    return m_values


### SIEMPRE HAY M ANTES DE N entonces busco los bloques antes de N y ahi me quedo

def define_block(cigar,blocks,cigartuples):
    bloques = []
    for x in [index for index,value in enumerate(cigar) if value =="N"]:
        m_index = cigar[:x].count("M") -1 
        bloques.append(blocks[m_index][1])
    return bloques

def define_block_two_N(cigar,blocks,m_values_index,cigartuples):
    bloques = []
    count = -1
    for x in [index for index,value in enumerate(cigar) if value =="N"]:
        count += 1
        if m_values_index == count:
            m_index = cigar[:x].count("M") -1 
            bloques.append(blocks[m_index][1])
        else:
            pass
    return bloques[0]

def closest(list, Number):
    """
    Para obtener el numero mas cercano a 500 y asi ver cual es la insercion mas probable
    """
    aux = []
    for valor in list:
        aux.append(abs(Number-valor))
    return aux.index(min(aux))    


def define_unique_maps(samfilete,te_list,df_te_info):
    TE_lost_bad_mapped = []
    TE_lost_bad_mapped_reasons = []
    header = ["chromosome","start","end","ID","frame","strand","te_number","te_family","te_len","predicted_len","te_chr_in_strain","te_start_in_strain","te_end_in_strain","cigar_S_bases","INREF","low_quality"]
    list_te = []
    for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
        if read.query_name in te_list:
            te_len = df_te_info[df_te_info["teid"]==read.query_name].te_len.values[0] # obtengo el largo del TE en cada genoma
            te_number = df_te_info[df_te_info["teid"]==read.query_name].te_number.values[0] # es el numero identificatorio del te en la cepa
            te_family = df_te_info[df_te_info["teid"]==read.query_name].te_family.values[0]
            pb_right2map = df_te_info[df_te_info["teid"]==read.query_name].pb_right2map.values[0]
            pb_left2map = df_te_info[df_te_info["teid"]==read.query_name].pb_left2map.values[0]
            te_start_in_strain  = df_te_info[df_te_info["teid"]==read.query_name].te_start.values[0]
            te_end_in_strain = df_te_info[df_te_info["teid"]==read.query_name].te_end.values[0]
            te_chr_in_strain = df_te_info[df_te_info["teid"]==read.query_name].seqname.values[0]
            cigar = read.cigarstring # guardo el CIGAR
            cigarN = cigar.count("N")
            cigar_S_bases = read.get_cigar_stats()[0][4] # cuantas bases son tipo softclipping
            low_quality = "N"
            if cigarN > 2:
                TE_lost_bad_mapped_reasons.append("N+2")
                TE_lost_bad_mapped.append(read.query_name)
            else:
                if read.is_reverse: # establezco el strand
                    strand = "-"
                else:
                    strand = "+"
                length = read.reference_length # Me quedo con el largo de la referencia
                mappingquality = read.mapping_quality # Me quedo con la calidad de mapeo
                cigarN_bases= read.get_cigar_stats()[0][3] # cuantas bases son tipo intron
                cigarM_bases = read.get_cigar_stats()[0][0] # cuantas bases son tipo match/missmatch
                blocks = read.get_blocks()
                n_values = [x[1] for index,x in enumerate(read.cigartuples) if x[0] == 3]
                m_values = get_m_values(read.cigartuples)
                n_index = [index for index,x in enumerate(read.cigartuples) if x[0] == 3]
                if cigarN == 1:
                    n_values =n_values[0] # hay un solo N hay un solo valor
                    te_reference_start = define_block(cigar,blocks,read.cigartuples)[0]
                    te_reference_end = define_block(cigar,blocks,read.cigartuples)[0] + n_values  # hay una sola coordenada
                    predicted_len = te_reference_end -te_reference_start # Seria el largo en la referencia
                    diff_len = abs(predicted_len - te_len)
                    if diff_len > 200:
                        if (m_values[0] > 510 or m_values[0] < 490) and (m_values[1] > 510 or m_values[1] < 490):
                            low_quality = "Could be a different insertion in the reference."
                        else:
                            low_quality = "Big difference between predicted size and te size"
                elif cigarN > 1:
                    if cigar_S_bases > 100:
                        TE_lost_bad_mapped_reasons.append("N+1 soft-clipping")
                        TE_lost_bad_mapped.append(read.query_name)
                    else:
                        m_values_index = m_values[:1] + m_values[-1:] # me quedo con las coordenadas de los extremos ya que el TE tiene que ser una de ellas ----TE--------TE- en este caso el extremo derecho es muy chico por lo que la region con un match mas cercano a 500 es mas creible
                        m_values_index = closest(m_values_index, 500)
                        te_reference_start = define_block_two_N(cigar,blocks,m_values_index,read.cigartuples)
                        te_reference_end = define_block_two_N(cigar,blocks,m_values_index,read.cigartuples) + n_values[m_values_index]
                        predicted_len = te_reference_end -te_reference_start # Seria el largo en la referencia
                        diff_len = abs(predicted_len - te_len)
                        if m_values_index == 1:
                            m_values_index = 2 # m_values_index+1 porque es la ultima coordenada la que me interesa y anteriormente las habia redifinido
                        if (m_values[m_values_index] > 510 or m_values[m_values_index] < 490): 
                            TE_lost_bad_mapped.append(read.query_name)
                            TE_lost_bad_mapped_reasons.append("N+1 location")
                            # las elimino porque son dificiles de definir...
                        else:
                            if diff_len >= 3000:
                                TE_lost_bad_mapped.append(read.query_name)
                                TE_lost_bad_mapped_reasons.append("N+1 size")
                                # la diferencia de tamaño es muy grande
                            else:
                                low_quality = "There is a second insertion in the reference"
                else:
                    if cigar_S_bases > 200:
                        TE_lost_bad_mapped.append(read.query_name)
                        TE_lost_bad_mapped_reasons.append("lost by soft_clipping")
                    else:
                        align_region = abs(read.reference_start -read.reference_end)
                        if align_region >= 500 and align_region <= 1200:
                            breakpoint = int(np.median([read.reference_start,read.reference_end]))
                            te_reference_start = breakpoint # defino un unico punto que es el breakpoint
                            te_reference_end = breakpoint
                            predicted_len = 0
                        else:
                            TE_lost_bad_mapped_reasons.append("N=0, not sure region")
                            TE_lost_bad_mapped.append(read.query_name)
                if read.query_name not in TE_lost_bad_mapped:
                    if te_chr_in_strain != read.reference_name:
                        low_quality = "different chromosome"  # la determinacion dio en otro cromosoma
                    if te_reference_start == te_reference_end:
                        INREF = "N"
                    else:
                        INREF = "Y"
                    list_te.append([read.reference_name,te_reference_start,te_reference_end,read.query_name,".",strand,te_number,te_family,te_len,predicted_len,te_chr_in_strain,te_start_in_strain,te_end_in_strain,cigar_S_bases,INREF,low_quality])
    samfilete.close()
    df_unique = pd.DataFrame(list_te,columns=header)
    return df_unique,TE_lost_bad_mapped,TE_lost_bad_mapped_reasons


def runcml(cml):
    """Ejecuta un comando"""
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    print(stderr)
    stdout = stdout.decode("utf-8")
    return stdout

def runminimap2(strain_TE_file,reference_genome_file,outfile,mode):
    """
    strain_TE_file: son secuencias de TE (o flanqueantes) del strain
    reference_genome_file: es el genoma de ISO1
    mode: complete o flanking (archivo fasta de regiones flanqueantes a TE o de los propios TEs)
    """
    if mode == "complete":
        commandLine = "minimap2 -ax asm5 -t 8  -uf -C5 " + reference_genome_file + " " + strain_TE_file
    elif mode == "flanking": 
        commandLine = "minimap2 -ax splice -G100k -B2 -t 8 -uf -C5  " + reference_genome_file + " " + strain_TE_file
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def bedtoolsclosest(reference_gene_bed,te_bed,outfile,mode):
    """
    reference_gene_bed: es el archivo de genes de la referencia preparado para cada genoma
    te_bed son las coordenadas de los TE de cada genoma en el genoma de referencia
    """
    if mode == "overlap":
        commandLine = "bedtools  closest -D ref -t all -a " + te_bed + " -b " + reference_gene_bed
    elif mode == "upstream": # -io ignoro overlap, ignoro -id downstream
        commandLine = "bedtools  closest -D ref -t all -io -id -a " + te_bed + " -b " + reference_gene_bed
    else:
        # -io ignoro overlap, ignoro -iu upstream
        commandLine = "bedtools  closest -D ref -t all -io -iu -a " + te_bed + " -b " + reference_gene_bed
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def parse_closest_files(closest_file,mode):
    header = ["chr","start","end","TEID","frame","strand","chr_1","start_1","end_1","gene","frame_1","strand_1","distance"]
    used_cols = ["chr","start","end","TEID","frame","strand","gene","distance"]
    df_bed = pd.read_table(closest_file,sep="\t",header=None,names=header)
    df_bed = df_bed[used_cols]
    if sum(df_bed["frame"] == ".") == 0:
       df_bed["frame"] = "."
    if mode == "overlap":
        df_bed = df_bed[df_bed["distance"]==0]
    return df_bed    


def intersect(strain_te,reference_te,outfile):
    """
    El objetivo de este paso es crear un archivo de interseccion entre los TEs del genoma del strain pasados al genoma de referencia, para poder determinar coincidencias
    """
    commandLine = "bedtools  intersect -a " + strain_te + " -b " + reference_te + " -wao "
    print(commandLine)
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def parse_intersect_files(intersection_file):
    header = ["chr","start","end","te_strain","frame","strand","chr_1","start_1","end_1","te_ref","frame_1","strand_1","intersect"]
    df_bed = pd.read_table(intersection_file,sep="\t",header=None,names=header)
    used_cols = ['te_strain',"te_ref",'intersect']
    df_bed = df_bed[used_cols]
    return df_bed


def intersect_euchromatin(strain_te_in_ref,reference_euchromatin,outfile):
    """
    El objetivo de este paso es crear un archivo de interseccion entre los TEs del genoma del strain pasados al genoma de referencia, para poder determinar coincidencias
    """
    commandLine = "bedtools  intersect -a " + strain_te_in_ref + " -b " + reference_euchromatin + " -wao "
    print(commandLine)
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile

def parse_intersect_files_euchromatin(intersection_file):
    header = ["chr","start","end","te_strain","frame","strand","chr_ref","start_ref","end_ref","intersect"]
    df_bed = pd.read_table(intersection_file,sep="\t",header=None,names=header)
    used_cols = ['te_strain','start_ref']
    df_bed = df_bed[used_cols]
    return df_bed



#### Analisis de TEs

## Cargo el log file 



logfile = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2]+ ".log"

### Primero que todo creo un directorio 

fo = open(logfile,"a")
log_line = "--- Starting TE Transfer from "+ sys.argv[2] + " to ISO1 ---\n"
fo.write(log_line)
fo.close()


### Primero veo si el mapeo de la secuencia del TE fue o no unico

### Veo si la secuencia completa de los TEs mapearon en un lugar unico




reference_genome_file = sys.argv[3]

te_complete_sequence =  sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_sequence.fasta"

te_complete_samfile = sys.argv[1] + "/" +  sys.argv[2] + "/TE_sequence_" + sys.argv[2] + "_2_ISO.sam"


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for complete TE sequences ---\n")
fo.close()

runminimap2(te_complete_sequence,reference_genome_file,te_complete_samfile,"complete")


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for complete TE sequences: Done ---\n")
fo.close()


samfilete = pysam.AlignmentFile(te_complete_samfile, "r")

list_te_complete = []

for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    cigar_S_bases = read.get_cigar_stats()[0][4] # cuantas bases son tipo softclipping
    list_te_complete.append([read.query_name,read.reference_name,read.reference_start,read.reference_end,read.is_unmapped,cigar_S_bases])

samfilete.close()


# Determino mapeos unicos

df_complete_seq_mapping = pd.DataFrame(list_te_complete,columns=["id","chr","start","end","unmapped","cigar_S_bases"])

te_complete_lunique = []
for name,grp in df_complete_seq_mapping.groupby("id"):
    if grp.unmapped.values[0]:
        pass
    elif len(grp) == 1:
        te_complete_lunique.append(grp.id.values[0])
    else:
        pass


te_complete_unique = len(te_complete_lunique) # 634

fo = open(logfile,"a")
fo.write("--- Determining the number of unique assingments for TEs complete sequence ---\n")
fo.close()


fo = open(logfile,"a")
log_line = "--- There are in total: " + str(te_complete_unique) + " ---\n"
fo.write(log_line)
fo.close()


# Luego cargo el resultado de minimap2 de las regiones flanqueantes a los TEs +- 500 (un archivo SAM) a un objeto pysam, para poder iterar cada read del resultado

## Defino todos los reads existentes y los llevo a un dataframe

te_flanking_samfile = sys.argv[1] + "/" + sys.argv[2] + "/TE_flanking_regions_" + sys.argv[2] + "_2_ISO.sam"

te_flanking_sequence = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_flanking_regions.fasta"


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for flanking TE sequences ---\n")
fo.close()

runminimap2(te_flanking_sequence,reference_genome_file,te_flanking_samfile,"flanking")


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for flanking TE sequences: Done ---\n")
fo.close()

samfilete = pysam.AlignmentFile(te_flanking_samfile, "r")

list_te = []
for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    list_te.append([read.query_name,read.reference_start,read.reference_end,read.is_unmapped])

samfilete.close()


# Determino mapeos unicos y no unicos, tambien no mapeados

df_to_define_type_of_mapping = pd.DataFrame(list_te,columns=["teid","start","end","unmapped"])

## Cargo la informacion de los TE (como fue su transferencia)

df_te_info = pd.read_table("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Results/"+ sys.argv[1] + "/" + sys.argv[1] + "_TE_info.tsv",sep="\t",header=0)

df_te_info["te_number"] = df_te_info["teid"].str.split("_").str[1:2].str[0] # dicatamino numero de te
df_te_info["te_family"] = ["_".join(ID.split("_")[2:]) for ID in df_te_info.teid.tolist()]

# df_te_info["teid"].str.split("_").str[2:].str[0]  # dicatamino familia



fo = open(logfile,"a")
log_line = "--- Total TEs annotated: "+ str(len(df_te_info))+" ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Total TEs removed prior transfer (both boundaries incomplete (neither reach 500 pb)): "+ str(len(df_te_info[df_te_info["2_transfer"]=="N"]))+" ---\n"
fo.write(log_line)
fo.close()


## Los unicos seran aquellos que mapean en un solo lugar y parten de secuencias con regiones flanqueantes completas

df_to_define_type_of_mapping = df_to_define_type_of_mapping.merge(df_te_info,on="teid")
df_to_define_type_of_mapping["end"] = df_to_define_type_of_mapping["end"].fillna(0).astype(int)

fo = open(logfile,"a")
log_line = "--- Determining the number of unique mapping, unique maping but incomplete match (TE has N or N in on of the flanking regions), multimapping or not mapped flanking regions ---\n"
fo.write(log_line)
fo.close()

lnotmapped = []
lnotunique = []
lunique = []
lunique_problems = []
for name,grp in df_to_define_type_of_mapping.groupby("teid"):
    if grp.start.values[0] == -1:
        lnotmapped.append(grp.teid.values[0])
    else:
        if len(grp) > 1:
            lnotunique.append(grp.teid.values[0])
        else:
            if grp.contig_boundaries.values[0] == "Y" or grp.has_n.values[0] == "Y":
                lunique_problems.append(grp.teid.values[0])
            else:
                lunique.append(grp.teid.values[0])
    

fo = open(logfile,"a")
log_line = "--- Number of unique mapping: " + str(len(lunique)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Number of unique maping but incomplete match: " + str(len(lunique_problems)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Number of multimapping: " + str(len(lnotunique)) + " ---\n"
fo.write(log_line)
fo.close()


fo = open(logfile,"a")
log_line = "--- Number of not mapped: " + str(len(lnotmapped)) + " ---\n"
fo.write(log_line)
fo.close()



# Vuelvo a iterar los read para definir donde mapean aquellas regiones que fueron definidas como unicas 

fo = open(logfile,"a")
log_line = "--- Defining unique mapping location ---\n"
fo.write(log_line)
fo.close()

samfilete = pysam.AlignmentFile(te_flanking_samfile, "r")

df_unique,TE_lost_bad_mapped,TE_lost_bad_mapped_reasons = define_unique_maps(samfilete,lunique,df_te_info)





fo = open(logfile,"a")
log_line = "--- Number of unique mapping transferred: " + str(len(df_unique)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Number of unique mapping not transferred: " + str(len(TE_lost_bad_mapped)) + " ---\n"
fo.write(log_line)
fo.close()

# Hago lo mismos con los que no pude conseguir la region completa pero mapearon de forma unica

fo = open(logfile,"a")
log_line = "--- Defining unique incomplete flanking regions location ---\n"
fo.write(log_line)
fo.close()

samfilete = pysam.AlignmentFile(te_flanking_samfile, "r")

df_unique_n_contig,TE_lost_bad_mapped_with_problems,TE_lost_bad_mapped_reasons_with_problems = define_unique_maps(samfilete,lunique_problems,df_te_info)


fo = open(logfile,"a")
log_line = "--- Number of unique incomplete flanking regions transferred: " + str(len(df_unique_n_contig)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Number of unique incomplete flanking regions not transferred: " + str(len(TE_lost_bad_mapped_with_problems)) + " ---\n"
fo.write(log_line)
fo.close()


TE_lost_bad_mapped += TE_lost_bad_mapped_with_problems
TE_lost_bad_mapped_reasons += TE_lost_bad_mapped_reasons_with_problems


### junto la info de las dos determinaciones

df_unique = pd.concat([df_unique,df_unique_n_contig])


fo = open(logfile,"a")
log_line = "--- Number of unique and unique-incomplete transferred: " + str(len(df_unique)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Number of unique and unique-incomplete not transferred: " + str(len(TE_lost_bad_mapped)) + " ---\n"
fo.write(log_line)
fo.close()

###############################################################################################
###############################################################################################

### Tengo que trabajar con los que mapearon con problemas


## Una vez definido los mapeos unicos trabajo con la lista de TEs donde obtengo MM o problemas

# Sumo a los no unicos lo mal mapeados y los no mapeados

fo = open(logfile,"a")
log_line = "--- Processing multimapping and not transferred regions in order to perform reallocation ---\n"
fo.write(log_line)
fo.close()


te_to_remap = lnotunique + TE_lost_bad_mapped
te_to_remap.sort() 

fo = open(logfile,"a")
log_line = "--- There are: " + str(len(te_to_remap)) +" TE's ---\n"
fo.write(log_line)
fo.close()

## Obtengo informacion de los read que entran en la categoria remap

samfilete = pysam.AlignmentFile(te_flanking_samfile, "r")

rows = []

for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    if read.query_name in te_to_remap:
        te_len = df_te_info[df_te_info["teid"]==read.query_name].te_len.values[0] # obtengo el largo del TE en cada genoma
        te_name = read.query_name
        pb_right2map = df_te_info[df_te_info["teid"]==read.query_name].pb_right2map.values[0]
        pb_left2map = df_te_info[df_te_info["teid"]==read.query_name].pb_left2map.values[0]
        te_start_in_strain  = df_te_info[df_te_info["teid"]==read.query_name].te_start.values[0]
        te_end_in_strain = df_te_info[df_te_info["teid"]==read.query_name].te_end.values[0]
        te_chr_in_strain = df_te_info[df_te_info["teid"]==read.query_name].seqname.values[0]
        cigar = read.cigarstring # guardo el CIGAR
        cigar_S_bases = read.get_cigar_stats()[0][4] # cuantas bases son tipo softclipping
        length = read.reference_length # Me quedo con el largo de la referencia
        mappingquality = read.mapping_quality # Me quedo con la calidad de mapeo
        cigarN_bases= read.get_cigar_stats()[0][3] # cuantas bases son tipo intron
        cigarM_bases = read.get_cigar_stats()[0][0] # cuantas bases son tipo match/missmatch
        blocks = read.get_blocks()
        n_values = [x[1] for index,x in enumerate(read.cigartuples) if x[0] == 3]
        m_values = get_m_values(read.cigartuples)
        n_index = [index for index,x in enumerate(read.cigartuples) if x[0] == 3]
        te_start = df_te_info[df_te_info["teid"] ==te_name].te_start.values[0]
        te_end = df_te_info[df_te_info["teid"] ==te_name].te_end.values[0]
        te_chr = df_te_info[df_te_info["teid"] ==te_name].seqname.values[0]
        if read.is_reverse: # establezco el strand
            strand = "-"
        else:
            strand = "+"
        diff_te_pos_inic = abs(te_start - read.reference_start)
        if length > 200 and te_chr_in_strain == read.reference_name:
                rows.append([te_name,te_len,read.reference_name,strand,read.reference_start,read.reference_end,te_chr_in_strain,te_start,te_end,diff_te_pos_inic,cigar,length,cigarM_bases,mappingquality,cigarN_bases,read.cigartuples,blocks,cigar_S_bases,read.flag])



df_to_refine_mm = pd.DataFrame(rows,columns=["te_name","te_len","reference_name","strand","reference_start","reference_end","te_chr_in_strain","te_start","te_end","diff_te_pos_inic","cigar","length","cigarM_bases","mappingquality","cigarN_bases","cigartuples","bloque","cigar_S_bases","flag"])

samfilete.close()

# Intento reducir los mapeos multiples utilizando la info de el mapeo completo y unico de secuencia, la mas cercana sera la original


rows=[]
for name,grp in df_to_refine_mm.groupby("te_name"):
    flag = True
    if len(grp) > 1:
        if name in te_complete_lunique:
            temp = df_complete_seq_mapping[df_complete_seq_mapping["id"]==name]
            if temp.chr.values[0] == grp.reference_name.values[0] and len(temp) ==1:
                most_near_pos = 1000000
                for index,row in grp.iterrows():
                    most_near_pos_i = abs(row.reference_start-temp.start.values[0])
                    if most_near_pos_i < most_near_pos:
                        most_near_pos = most_near_pos_i
                        real_pos = row.reference_start
                if most_near_pos != 1000000:
                    rows.append(grp[grp["reference_start"]==real_pos].values[0].tolist())
                    flag = False
    if flag:
        for index,row in grp.iterrows():
            rows.append(row.tolist())

df_to_refine_mm = pd.DataFrame(rows,columns=df_to_refine_mm.columns.tolist())


te_set = set(df_to_refine_mm.te_name.tolist())

samfilete = pysam.AlignmentFile(te_flanking_samfile, "r")
te_mmap_count = defaultdict(int)

for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    if read.query_name in te_set:
        te_mmap_count[read.query_name] += 1

samfilete.close()



### En caso que existan mapeos suplementarios me quedo con los que mapearon con mas bases M y si el numero es el mismo me quedo con los casos que la pos de inicio del te es mas cercano a la ref



df_to_refine_mm = df_to_refine_mm.sort_values(["te_name","cigarM_bases","diff_te_pos_inic"],ascending=[True,False,True])

df_to_refine_mm = df_to_refine_mm.drop_duplicates("te_name",keep="first")
df_to_refine_mm = df_to_refine_mm.sort_values(["reference_name","reference_start","reference_end"],ascending=[True,True,True])
df_to_refine_mm = df_to_refine_mm.reset_index(drop=True)


TE_not_assigned = []
TE_not_assigned_reason = []

header = ["chromosome","start","end","ID","frame","strand","te_number","te_family","te_len","predicted_len","te_chr_in_strain","te_start_in_strain","te_end_in_strain","cigar_S_bases","INREF","low_quality"]
list_te = []
for index,row in df_to_refine_mm.iterrows():
    name = row.te_name 
    if row.diff_te_pos_inic > row.reference_start:
        TE_not_assigned.append(name)
        TE_not_assigned_reason.append("diff_start_pos")
    else:
        INREF = "Undefined"
        predicted_len = 0
        cigar = row.cigar
        cigarN = cigar.count("N")
        low_quality = "N"
        te_len = df_te_info[df_te_info["teid"]==name].te_len.values[0] # obtengo el largo del TE en cada genoma
        te_number = df_te_info[df_te_info["teid"]==name].te_number.values[0]
        te_family = df_te_info[df_te_info["teid"]==name].te_family.values[0]
        te_chr_in_strain  = df_te_info[df_te_info["teid"]==name].seqname.values[0]
        te_start_in_strain  = df_te_info[df_te_info["teid"]==name].te_start.values[0]
        te_end_in_strain = df_te_info[df_te_info["teid"]==name].te_end.values[0]
        te_chr_in_strain = df_te_info[df_te_info["teid"]==name].seqname.values[0]
        chromosome = row.reference_name
        cigar_S_bases = row.cigar_S_bases
        cigarN_bases = row.cigarN_bases
        cigartuples = row.cigartuples
        length = row.length
        strand = row.strand
        reference_start = row.reference_start
        reference_end = row.reference_end
        blocks = row.bloque
        n_values = [x[1] for index,x in enumerate(cigartuples) if x[0] == 3]
        m_values = get_m_values(cigartuples)
        n_index = [index for index,x in enumerate(cigartuples) if x[0] == 3]
        if cigarN == 1:
            INREF = "Y"
            n_values =n_values[0] # hay un solo N hay un solo valor
            te_reference_start = define_block(cigar,blocks,cigartuples)[0]
            te_reference_end = define_block(cigar,blocks,cigartuples)[0] + n_values  # hay una sola coordenada
            predicted_len = te_reference_end -te_reference_start # Seria el largo en la referencia
            diff_len = abs(predicted_len - te_len)
            percent_size = int(te_len*100/predicted_len)
            if percent_size >= 95 and percent_size <= 105:
                if diff_len > 200:
                    INREF = "Undefined"
                    low_quality = "Big difference between predicted size and te size"
            else:
                TE_not_assigned.append(name)
                TE_not_assigned_reason.append("Multimapping region, size dont match. porbably nested TE or in tandem")
        elif cigarN > 1:
            INREF = "Y"
            m_values_index = m_values[:1] + m_values[-1:] # me quedo con las coordenadas de los extremos ya que el TE tiene que ser una de ellas ----TE--------TE- en este caso el extremo derecho es muy chico por lo que la region con un match mas cercano a 500 es mas creible
            m_values_index = closest(m_values_index, 500)
            te_reference_start = define_block_two_N(cigar,blocks,m_values_index,cigartuples)
            te_reference_end = define_block_two_N(cigar,blocks,m_values_index,cigartuples) + n_values[m_values_index]
            predicted_len = te_reference_end -te_reference_start # Seria el largo en la referencia
            diff_len = abs(predicted_len - te_len)
            percent_size = int(te_len*100/predicted_len)
            if percent_size >= 95 and percent_size <= 105:
                if m_values_index == 1:
                    m_values_index = 2 # m_values_index+1 porque es la ultima coordenada la que me interesa y anteriormente las habia redifinido
                if (m_values[m_values_index] > 510 or m_values[m_values_index] < 490): 
                    TE_not_assigned.append(name)
                    TE_not_assigned_reason.append("Multimapping region and N+1 location. Probably nested TE or in tandem")
                    # las elimino porque son dificiles de definir...
                else:
                    low_quality = "Multimapping region"
            else:
                TE_not_assigned.append(name)
                TE_not_assigned_reason.append("Multimapping region and N+1 location. Probably nested TE or in tandem")
        else:
            INREF = "Undefined"
            align_region = abs(reference_start - reference_end)
            if strand == "+":
                te_reference_end = reference_end + te_len + 500
                te_reference_start = reference_end
            else:
                te_reference_start = reference_end - 100
                te_reference_end = reference_end + 100
            if align_region >= 250:
                if cigar_S_bases > 450:
                    low_quality = "Incomplete region match"
                else:  
                    low_quality = "Multimapping region"
            else:
                TE_not_assigned.append(name)
                TE_not_assigned_reason.append("Multimapping region can not transfer")
        if name not in TE_not_assigned:
            if chromosome != te_chr_in_strain:
                if low_quality == "N":
                    low_quality = "different chromosome"
                else:
                    low_quality = "different chromosome. " + low_quality
            list_te.append([chromosome,te_reference_start,te_reference_end,name,".",strand,te_number,te_family,te_len,predicted_len,te_chr_in_strain,te_start_in_strain,te_end_in_strain,cigar_S_bases,INREF,low_quality])    


df_mm_2_unique = pd.DataFrame(list_te,columns=header)


df_mm_2_unique = df_mm_2_unique.sort_values(["chromosome","start","end"],ascending=[True,True,True])
df_mm_2_unique = df_mm_2_unique.reset_index(drop=True)


# for index,row in df_mm_2_unique.iterrows():
#    print(row.ID,row.chromosome,row.start)

fo = open(logfile,"a")
log_line = "--- TE transferred in the reallocation: " + str(len(df_mm_2_unique)) + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- TE not transferred in the reallocation: " + str(len(TE_not_assigned)) + " ---\n"
fo.write(log_line)
fo.close()


### ahora concateno la informacion de los unicos con los mm

df_te = pd.concat([df_unique,df_mm_2_unique])

df_te = df_te.sort_values(["chromosome","start","end"])

df_te = df_te.reset_index(drop=True)


fo = open(logfile,"a")
log_line = "--- Total TE transferred: " + str(len(df_te)) + " ---\n"
fo.write(log_line)
fo.close()


# Voy a intentar obtener informacion de los que no mapearon
# para eso voy a buscar informacion 


### Agrego info sobre nested y tandem de cada TE


df_nested_tandem = pd.read_table( sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_nested_tandem_info.tsv",sep="\t",header=0)

te_transferred = set(df_te.ID.tolist()) # Veo cuales fueron los TE transferidos - 1953
list_te_to_redefine = lnotmapped + TE_not_assigned # Veo cuales fueron los TE transferidos

df_nested = df_nested_tandem[df_nested_tandem["inside_of"]!="."][["ID","inside_of"]]

nested_not_predicted_but_parent_found = [te for te in set(df_nested.ID.tolist()) if te not in te_transferred and df_nested[df_nested["ID"] == te].inside_of.values[0] in  te_transferred]

set_nested = set(nested_not_predicted_but_parent_found) # sumo el conjunto de nested que quiero mantener


fo = open(logfile,"a")
log_line = "--- TE nested with TE parent already transferred: " + str(len(set_nested)) + " ---\n"
fo.write(log_line)
fo.close()

assinged_nested_te = []

for te in set_nested:
    tandem_up_in_ref = ""
    tandem_down_in_ref = ""
    te_chr_in_strain = df_te_info[df_te_info["teid"]==te].seqname.values[0] # cromosoma
    te_start_in_strain = df_te_info[df_te_info["teid"]==te].te_start.values[0] # inicio en strain
    te_end_in_strain = df_te_info[df_te_info["teid"]==te].te_end.values[0] # fin en strain
    te_family = te.split("_")[2:][0] # recupero familia
    te_number = te.split("_")[1] # recupero numero identificatorio
    frame = "." # defino el frame
    strand = df_te_info[df_te_info["teid"]==te].te_strand.values[0]
    te_len = df_te_info[df_te_info["teid"]==te].te_len.values[0]
    cigar_S_bases = "NA"
    inside_of= df_nested_tandem[df_nested_tandem["ID"]==te].inside_of.values[0]
    # utilizaremos dos formas distintas de resolverlo, si inside of esta o no en la refrencia y el nested te no
    # . fijaremos la posicion a partir de la distancia desde el incio del nested al inicio de parent. Fijaremos como un breakpoint
    inside_of_inref = df_te[df_te["ID"]==inside_of].INREF.values[0]
    if inside_of_inref == "N": ## expando en 5 bases el parent y fijo el breakpoint para el nested
        INREF = inside_of_inref
        predicted_len = 0
        breakpoint = df_te[df_te["ID"]==inside_of].start.values[0]
        chromosome = df_te[df_te["ID"]==inside_of].chromosome.values[0]
        parent_te = df_te[df_te["ID"]==inside_of].values[0].tolist()
        parent_te[1] = parent_te[1] - 5 # EXPANDO EN 5 EL INICIO Y FIN DEL PADRE
        parent_te[2] = parent_te[2] + 5 # EXPANDO EN 5 EL INICIO Y FIN DEL PADRE
        assinged_nested_te.append(parent_te)
        low_quality = "Nested TE assigned"
        assinged_nested_te.append([chromosome,breakpoint,breakpoint,te,frame,strand,te_number,te_family,te_len,0,te_chr_in_strain,te_start_in_strain,te_end_in_strain,"NA",INREF,low_quality])
        

df_assign_nested = pd.DataFrame(assinged_nested_te,columns=df_te.columns.tolist())

if len(df_assign_nested) > 0:
    te_redefined = int(len(df_assign_nested)/2)
else:
    te_redefined = 0


fo = open(logfile,"a")
log_line = "--- TE nested transferred: " + str(te_redefined) + " ---\n"
fo.write(log_line)
fo.close()


### ahora concateno la informacion de los ya definidos con los nested

df_te = pd.concat([df_te,df_assign_nested])

df_te = df_te.sort_values(["chromosome","start","end"])

df_te = df_te.drop_duplicates("ID",keep="first")

df_te = df_te.reset_index(drop=True)

# for index,row in df_te.iterrows():
#    print(row.ID,row.chromosome,row.start,te_mmap_count[row.ID])



te_transferred = set(df_te.ID.tolist()) # Veo cuales fueron los TE transferidos


fo = open(logfile,"a")
log_line = "--- TE transferred: " + str(len(te_transferred)) + " ---\n"
fo.write(log_line)
fo.close()


## Me fijo cuantos no pase

te_not_transferred = [te for te in set(df_complete_seq_mapping.id) if te not in te_transferred]

fo = open(logfile,"a")
log_line = "--- TE not transferred: " + str(len(te_not_transferred)) + " ---\n"
fo.write(log_line)
fo.close()


### Agrego informacion sobre si la secuencia mapeo unica y si la region flanqueante lo hizo 

flanking_unique_mapping = []
for index,row in df_te.iterrows():
    if row.ID in lunique:
        flanking_unique_mapping += ["Y"]
    else:
        flanking_unique_mapping += ["N"]

complete_seq_mapping_unique = []
for index,row in df_te.iterrows():
    if row.ID in te_complete_lunique:
        if row.chromosome == df_complete_seq_mapping[df_complete_seq_mapping["id"]==row.ID].chr.values[0] and (df_complete_seq_mapping[df_complete_seq_mapping["id"]==row.ID].cigar_S_bases.values[0] < 100):
            complete_seq_mapping_unique += ["Y"]
        else:
            complete_seq_mapping_unique += ["N"]
    else:
        complete_seq_mapping_unique += ["N"]


## Agrego info de flanqueante y complete seq


df_te["flanking_unique_mapping"] = flanking_unique_mapping

     
df_te["sequence_unique_mapping"] = complete_seq_mapping_unique    


### Me fijo si hay tes con las mismas coordenadas para eso defino un id con chr start y end
 
df_te["group_id"] = ["-".join(map(str,x)) for x in df_te[["chromosome","start","end"]].values.tolist()]

count = 0
te_merged = []
te_merged_ids = []
for name,grp in df_te.groupby("group_id"):
    number_te_merged = 1 
    if len(grp) > 1:
        number_te_merged = len(grp)
        grp = grp.sort_values(["low_quality"])
        te_merged_ids.append("&".join(grp.ID.values.tolist()))
    te_merged.append(grp.values[0].tolist() + [number_te_merged])

df_te = pd.DataFrame(te_merged,columns=df_te.columns.tolist() + ["number_te_merged"])


fo = open(logfile,"a")
log_line = "--- TE regions merged: " + str(len(te_merged_ids)) + " ---\n"
fo.write(log_line)
fo.close()

n_te_merged = sum([len(te.split("&")) for te in te_merged_ids])


fo = open(logfile,"a")
log_line = "--- TE merged: " + str(n_te_merged) + " ---\n"
fo.write(log_line)
fo.close()

te_merged_line = [te.split("&")[0] + "\t" + " - ".join(te.split("&")[1:]) for te in te_merged_ids]

te_merged_line = "\n".join(te_merged_line)

merged_tes_file = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_merged.tsv"

fo = open(merged_tes_file,"w")
fo.write(te_merged_line)
fo.close()

## Agrego info de la region del TE

df_te = df_te.merge(df_te_info[["teid","has_n","contig_boundaries"]],left_on="ID",right_on="teid")

df_te =df_te.drop("teid",axis=1)
df_te = df_te.drop("group_id",axis=1)

df_te = df_te.sort_values(["chromosome","start","end"])

df_te = df_te.reset_index(drop=True)


fo = open(logfile,"a")
log_line = "--- TE transferred without coordinate duplication: " + str(len(df_te)) + " ---\n"
fo.write(log_line)
fo.close()

df_te_backup = df_te

df_te_2_bed =df_te[['chromosome', 'start', 'end', 'ID', 'frame', 'strand']]

df_te_2_bed.to_csv(sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_transfer.bed",sep="\t",header=False,index=False)

###########################################################################################################################################################################################################################################################################################################


### Chequeo ubicacion con respecto a posicion de genes


# Voy a hacer 3 ejecuciones del closet, uno para sacar los que solapan con genes, otro para determinar los up y otro los down


# UTILIZO UNICAMENTE LOS GENES QUE NO RESULTARON PROBLEMATICOS



strain_gene_info = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer_info.tsv" 

df_gene_strain_info= pd.read_table(strain_gene_info,sep="\t",header=0)

softclipping = []
for index, row in df_gene_strain_info.iterrows():
    if "S" in row.cigar:
        softclipping += ["Y"]
    else:
        softclipping += ["N"]

df_gene_strain_info["softclipping"] = softclipping

genes_remove_softclipping = df_gene_strain_info[((df_gene_strain_info["percent_size"]<90) | (df_gene_strain_info["percent_size"]>110)) & (df_gene_strain_info["softclipping"] == "Y") & (df_gene_strain_info["diff_ref"]>=500)].gene.tolist()



strain_gene_bed = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer.bed" # Anotacion de los genes del strain

df_gene_strain= pd.read_table(strain_gene_bed,sep="\t",header=None,names=["chr","start","end","gene","frame","strand","MM","I","C","S"])

genes_good_for_transfer = []
for index,row in df_gene_strain.iterrows():
    if row.MM == "NMM" and row.I == "NI" and row.C == "NC" and row.S == "N":
        if row.gene not in genes_remove_softclipping:
            genes_good_for_transfer.append(row.gene)

len(genes_good_for_transfer) # 13259

df_gene_strain = df_gene_strain[df_gene_strain["gene"].isin(genes_good_for_transfer)][["chr","start","end","gene","frame","strand"]]

df_gene_strain = df_gene_strain.sort_values(["chr","start","end"])

df_gene_strain = df_gene_strain.reset_index(drop=True)

df_gene_strain.to_csv(sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer_for_comparision.bed",sep="\t",header=False,index=False)


reference_gene_bed = sys.argv[1] + "/"  + sys.argv[2] + "/Reference_gene_for_" + sys.argv[2] + "_transfer.bed"  # estas son los genes de REF que estan presentes en la anotacion de la cepa

df_gene_ref= pd.read_table(reference_gene_bed,sep="\t",header=None,names=["chr","start","end","gene","frame","strand"])


df_gene_ref = df_gene_ref[df_gene_ref["gene"].isin(genes_good_for_transfer)][["chr","start","end","gene","frame","strand"]]

df_gene_ref = df_gene_ref.sort_values(["chr","start","end"])

df_gene_ref = df_gene_ref.reset_index(drop=True)

df_gene_ref.to_csv(sys.argv[1] + "/"  + sys.argv[2] + "/Reference_gene_for_" + sys.argv[2] + "_transfer_for_comparision.bed",sep="\t",header=False,index=False)


fo = open(logfile,"a")
log_line = "--- Adding information about the conservation of the synteny for each TE. Comparing the genes up an down of a TE in the genome of the strain to the genes up an down of a TE in the reference genome once TE were transferred ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
log_line = "--- Only genes without tags (no inversion (or related), no multimapping (or related), no translocation to other chromosome (or related)). In total : " + str(len(df_gene_ref)) + " genes considered ---\n"
fo.write(log_line)
fo.close()


###################################################################################################
###################################################################################################

# Defino la interaccion de cada TE del genoma de interes en el genoma de referencia (utilizando las coordenadas de la transferencia)

reference_gene_bed = sys.argv[1] + "/"  + sys.argv[2] + "/Reference_gene_for_" + sys.argv[2] + "_transfer_for_comparision.bed"  # estas son los genes de REF que estan presentes en la anotacion de la cepa

te_bed_strain = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_transfer.bed" # estas son las coordenadas de los TEs de la cepa llevadas al genoma de referencia


overlap = bedtoolsclosest(reference_gene_bed,te_bed_strain,"overlap.txt","overlap")
upstream = bedtoolsclosest(reference_gene_bed,te_bed_strain,"upstream.txt","upstream")
downstream = bedtoolsclosest(reference_gene_bed,te_bed_strain,"downstream.txt","downstream")

### Utilizo la tabla df_te para recuperar la familia

te_bed_strain = pd.read_table(te_bed_strain,sep="\t",header=None)

te_bed_strain.columns = ['chromosome', 'start', 'end', 'ID', 'frame', 'strand']


# Creo un diccionario para guardar la info de los te (genes up, down y overlap)
te_dict_strain_ref = {te:[df_te[df_te["ID"]==te].te_family.values[0]]+[""]*6 for te in te_bed_strain.ID.unique().tolist()}


# te_dict_strain_ref  significa que utilice los te de strain llevados al genoma de referencia por lo que uso los genes del genoma de referencia

# Defino los overlap 

df_overlap = parse_closest_files(overlap,"overlap")

te_overlap = df_overlap.TEID.unique().tolist()

for name,grp in df_overlap.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = 0
    te_dict_strain_ref[name][3:5] = [gene,distance]

for te in te_dict_strain_ref.keys():
    if te not in te_overlap:
        te_dict_strain_ref[te][3:5] = ["-","-"]


# Defino los up 

df_up = parse_closest_files(upstream,"upstream")

te_up = df_up.TEID.unique().tolist()

for name,grp in df_up.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = abs(grp.distance.values[0])
    te_dict_strain_ref[name][1:3] = [gene,distance]


# Defino los down 

df_down = parse_closest_files(downstream,"downstream")

te_down = df_down.TEID.unique().tolist()

for name,grp in df_down.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = abs(grp.distance.values[0])
    te_dict_strain_ref[name][5:7] = [gene,distance]

          


###############################################################################################
###############################################################################################

# Analizo la interaccion de los TE con los genes en el strain

strain_gene_bed = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer_for_comparision.bed" # Anotacion de los genes del strain

te_original_strain_bed = sys.argv[4]


overlap = bedtoolsclosest(strain_gene_bed,te_original_strain_bed,"overlap.txt","overlap")
upstream = bedtoolsclosest(strain_gene_bed,te_original_strain_bed,"upstream.txt","upstream")
downstream = bedtoolsclosest(strain_gene_bed,te_original_strain_bed,"downstream.txt","downstream")

df_te_original_strain_bed = pd.read_table(te_original_strain_bed,sep="\t",header=None)

df_te_original_strain_bed.columns = ['chromosome', 'start', 'end', 'ID', 'len', 'strand']

### Aca no puedo Utilizar la tabla df_te para recuperar la familia porque no estan colapsados

# Creo un diccionario para guardar la info de los te (genes up, down y overlap)

te_dict_strain_strain = {te:["".join(te.split("_")[2:])]+[""]*6 for te in df_te_original_strain_bed.ID.unique().tolist()}

# te_dict_strain_strain  significa que utilice los te y los genes del strain en su propio genoma
# Defino los overlap 

df_overlap = parse_closest_files(overlap,"overlap")

te_overlap = df_overlap.TEID.unique().tolist()

for name,grp in df_overlap.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = 0
    te_dict_strain_strain[name][3:5] = [gene,distance]

for te in te_dict_strain_strain.keys():
    if te not in te_overlap:
        te_dict_strain_strain[te][3:5] = ["-","-"]


# Defino los up 

df_up = parse_closest_files(upstream,"upstream")

te_up = df_up.TEID.unique().tolist()

for name,grp in df_up.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = abs(grp.distance.values[0])
    te_dict_strain_strain[name][1:3] = [gene,distance]


# Defino los down 

df_down = parse_closest_files(downstream,"downstream")

te_down = df_down.TEID.unique().tolist()

for name,grp in df_down.groupby("TEID"):
    if len(grp)>1:
        gene = "-".join(grp.gene.values.tolist())
    else:
        gene = grp.gene.values[0]
    distance = abs(grp.distance.values[0])
    te_dict_strain_strain[name][5:7] = [gene,distance]


## ahora el diccionario tiene toda la informacion necesaria

### Chequeo ambos dicts para poner la informacion de la sintenia en el genoma original y en la transferencia


rows = []
for key,value in te_dict_strain_strain.items():
    if key in te_dict_strain_ref.keys():
        rows.append([key]+value+te_dict_strain_ref[key][1:])


header = ["ID","Family","upstream_gene_strain","upstream_dist_gene_strain","overlap_gene_strain","overlap_dist_gene_strain","downstream_gene_strain","downstream_dist_gene_strain","upstream_gene_ref","upstream_dist_gene_ref","overlap_gene_ref","overlap_dist_gene_ref","downstream_gene_ref","downstream_dist_gene_ref"]

df_te_gene_info = pd.DataFrame(rows,columns=header)

df_te_gene_info = df_te_gene_info[["ID","Family","upstream_gene_strain","upstream_dist_gene_strain","upstream_gene_ref","upstream_dist_gene_ref","overlap_gene_strain","overlap_dist_gene_strain","overlap_gene_ref","overlap_dist_gene_ref","downstream_gene_strain","downstream_dist_gene_strain","downstream_gene_ref","downstream_dist_gene_ref"]]

# Determino si la sintenia se mantiene o se pierde y como 
# cambian las distancia del TE a los genes en cada genoma

synteny_info = []
for index,row in df_te_gene_info.iterrows():
    flag =True
    inref = df_te[df_te["ID"] == row.ID].INREF.values[0]
    strain_reference_synteny = "MANTAINED"
    if row.upstream_gene_strain != row.upstream_gene_ref or row.downstream_gene_strain != row.downstream_gene_ref or row.overlap_gene_strain != row.overlap_gene_ref:
        strain_reference_synteny = "LOSS"
        difference_distance_up = "NA"
        difference_distance_down = "NA"
    else:
        difference_distance_up = abs(row.upstream_dist_gene_strain - row.upstream_dist_gene_ref)
        difference_distance_down = abs(row.downstream_dist_gene_strain - row.downstream_dist_gene_ref)
    synteny_info.append([row.ID,inref,strain_reference_synteny,difference_distance_up,difference_distance_down])
        

df_synteny_info = pd.DataFrame(synteny_info,columns=["ID","INREF","SYNTENY","DIST_UP_DIFF","DIST_DOWN_DIFF"]) 


fo = open(logfile,"a")
log_line = "--- Synteny information added ---\n"
fo.write(log_line)
fo.close()


###############################################################################################
###############################################################################################

## Ahora tengo que ver si los TE presentes en la referencia coinciden con los predichos en la transferencia

fo = open(logfile,"a")
log_line = "--- Intersection of TE's transferred to the reference with the TE's of the reference genome ---\n"
fo.write(log_line)
fo.close()


reference_te_bed = sys.argv[4] + "/ISO1_TE_FlyBase_Repet.bed"

reference_te_bed_info = sys.argv[4] + "/ISO1_TE_FlyBase_Repet_info.tsv"

df_ref_te_info = pd.read_table(reference_te_bed_info,sep="\t",header=0)

df_ref_te_info = df_ref_te_info[['DB', 'ID', 'ID_REPET', 'TE_LEN', 'FAMILY']]
df_ref_te_info.columns = ['DB', 'ID_REF', 'ID_REF_ALIAS', 'TE_LEN', 'FAMILY']



strain_2_reference_te_bed = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_transfer.bed" # estas son las coordenadas de los TEs de RAL-177 llevadas al genoma de referencia

intersect_te_strain_ref = intersect(strain_2_reference_te_bed,reference_te_bed,"intersect_te_strain_ref.txt")

df_intersect_te_strain_ref = parse_intersect_files(intersect_te_strain_ref)

rows = []
for index,row in df_intersect_te_strain_ref.iterrows():
    if row.te_ref == ".":
        row.intersect = -1
    rows.append(row)

df_intersect_te_strain_ref = pd.DataFrame(rows,columns=df_intersect_te_strain_ref.columns.tolist())

df_synteny_intersect_info = df_synteny_info.merge(df_intersect_te_strain_ref,left_on="ID",right_on="te_strain",how="left")

df_synteny_intersect_info = df_synteny_intersect_info.merge(df_ref_te_info,left_on="te_ref",right_on="ID_REF",how="left")


df_synteny_intersect_info = df_synteny_intersect_info.drop("ID_REF",axis=1)
df_synteny_intersect_info = df_synteny_intersect_info.drop("te_strain",axis=1)

df_synteny_intersect_info["TE_LEN"] = df_synteny_intersect_info["TE_LEN"].fillna(0).astype(int)
df_synteny_intersect_info["ID_REF_ALIAS"] = df_synteny_intersect_info["ID_REF_ALIAS"].fillna(".")
df_synteny_intersect_info["DB"] = df_synteny_intersect_info["DB"].fillna(".")
df_synteny_intersect_info["FAMILY"] = df_synteny_intersect_info["FAMILY"].fillna(".")


### Comparo y veo si coniciden las familias, como puede haber mas de una interseccion elimino los duplicados


backup = df_synteny_intersect_info
rows = []
for index,row in df_synteny_intersect_info.iterrows():
    comment = "NO"
    strain_te_fam = "_".join(map(str,row.ID.split("_")[2:]))
    #row.ID.split("_")[2:][0]
    if row.INREF == "N" and row.intersect >= 0:
        if strain_te_fam == row.FAMILY:
            row.INREF = "Y"
            comment = "Originally was not detected"
        else:
            comment = "WRONG FAMILY"
    elif row.INREF == "Y":
        if row.intersect ==-1:
            comment = "TE NOT ANNOTATED IN REF"
            row.INREF = "Undefined"
        elif strain_te_fam != row.FAMILY:
            comment = "WRONG FAMILY"
            row.INREF = "N"
    elif row.INREF == "Undefined":
        if row.intersect >= 0:
            if strain_te_fam == row.FAMILY:
                row.INREF = "Y"
            else:
                comment = "WRONG FAMILY"
                row.INREF = "N"
        else:
            row.INREF = "N"
    rows.append(row.tolist() + [comment])

df_synteny_intersect_info = pd.DataFrame(rows,columns=df_synteny_intersect_info.columns.tolist()+["comments"])




## Elimino duplicados porque pudieron caer en regiones donde hay mas de un TE en la misma zona
## y no coincide la familia


fo = open(logfile,"a")
log_line = "--- Eliminating association of transposable elements where there was more than one match and at leastt one of them maintains the family. We remove the wrong associations ---\n"
fo.write(log_line)
fo.close()


df_synteny_intersect_info = df_synteny_intersect_info.sort_values(["ID","comments"])


rows = []
for name,grp in df_synteny_intersect_info.groupby("ID"):
    if len(grp) > 1:
        if grp.comments.tolist().count("WRONG FAMILY") == len(grp):
            rows.append(grp.values[0].tolist()[:5] + ['.', 0, '.', '.', '.', 0, 'NO'])
        elif grp.comments.tolist().count("NO") >= 1:
            for index,row in grp[grp["comments"]=="NO"].iterrows():
                rows.append(row.tolist())
        else:
            rows.append(grp.values[0].tolist())
    else:
        rows.append(grp.values[0].tolist())


df_synteny_intersect_info = pd.DataFrame(rows,columns=df_synteny_intersect_info.columns.tolist())




df_te = df_te.merge(df_synteny_intersect_info,on="ID")

df_te = df_te[['chromosome', 'start', 'end', 'ID', 'frame', 'strand','number_te_merged','INREF_y', 'te_number','te_family', 'te_len', 'predicted_len','te_chr_in_strain','te_start_in_strain', 'te_end_in_strain', 'has_n', 'contig_boundaries','flanking_unique_mapping', 'sequence_unique_mapping', 'te_ref', 'intersect', 'DB','ID_REF_ALIAS', 'TE_LEN', 'FAMILY', 'SYNTENY','DIST_UP_DIFF', 'DIST_DOWN_DIFF','low_quality', 'comments']]


df_te.columns = ['chromosome', 'start', 'end', 'id', 'frame', 'strand','te_merged','inref', 'te_number_id','te_strain_family', 'te_strain_len', 'predicted_len_in_ref','te_chr_in_strain','te_start_in_strain', 'te_end_in_strain', 'te_has_n_in_strain', 'te_near_contig_boundaries_in_strain','te_flanking_unique_mapping_in_ref', 'sequence_unique_mapping', 'id_ref', 'intersection_te_strain_te_ref', 'DB','id_ref_alias', 'te_ref_len', 'te_ref_family', 'te_maintain_gene_synteny_in_transfer','distance_difference_up_gene', 'distance_difference_down_gene','problems_in_assignment', 'other_comments']


### Agrego info sobre nested y tandem de cada TE


df_nested_tandem = pd.read_table(sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_nested_tandem_info.tsv",sep="\t",header=0)

df_te = df_te.merge(df_nested_tandem,left_on="id",right_on="ID")

### agrego info sobre heterocromatina

te_strain_transferred_to_ref =sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_transfer.bed"

euchromatin_bed = sys.argv[4] + "/dmel-all-euchromatin-region.bed"


euchromatin_intersect = intersect_euchromatin(te_strain_transferred_to_ref,euchromatin_bed,"euchromatin_intersect.txt")

df_euchromatin = parse_intersect_files_euchromatin(euchromatin_intersect)

tes_heterochromatin = set(df_euchromatin[df_euchromatin["start_ref"] == -1].te_strain)

fo = open(logfile,"a")
log_line = "--- TE's transferred to heterochromatin regions: " + str(len(tes_heterochromatin)) + " ---\n"
fo.write(log_line)
fo.close()


rows=[]
for index,row in df_te.iterrows():
    if row.id in tes_heterochromatin:
        region = "heterochromatin"
    else:
        region = "euchromatin"
    rows.append(row.tolist() + [region])


df_te = pd.DataFrame(rows,columns=df_te.columns.tolist()+["region"])



df_te = df_te[['chromosome', 'start', 'end', 'id', 'te_number_id','frame', 'strand', 'te_merged',
        'te_strain_family', 'te_strain_len', 'inref','predicted_len_in_ref', 
       'id_ref', 'intersection_te_strain_te_ref', 'DB', 'id_ref_alias',
       'te_ref_len', 'te_ref_family', 'te_maintain_gene_synteny_in_transfer',
       'distance_difference_up_gene', 'distance_difference_down_gene',
       'problems_in_assignment', 'other_comments','te_chr_in_strain', 'te_start_in_strain',
       'te_end_in_strain', 'te_has_n_in_strain',
       'te_near_contig_boundaries_in_strain', 'inside_of',
       'partial_nested', 'tandem_up', 'tandem_down', 'te_flanking_unique_mapping_in_ref', 'sequence_unique_mapping','region']]


te_transferred_info =  sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_transferred_2_ISO_info.tsv"

df_te.to_csv(te_transferred_info,sep="\t",header=True,index=False)

fo = open(logfile,"a")
log_line = "--- TE transferred information in: " + te_transferred_info + " ---\n"
fo.write(log_line)
fo.close()

os.remove("upstream.txt")
os.remove("downstream.txt")
os.remove("overlap.txt")
os.remove("euchromatin_intersect.txt")
os.remove("intersect_te_strain_ref.txt")

fo = open(logfile,"a")
log_line = "--- " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + ": Done ---\n"
fo.write(log_line)
fo.close()


fo = open(logfile,"a")
log_line = "--- " + sys.argv[1] + ": Done ---\n"
fo.write(log_line)
fo.close()


############################################################################################################################################################################################################################################################################################################

