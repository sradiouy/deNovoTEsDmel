#!/usr/bin/env python -W ignore::DeprecationWarning
# encoding: utf-8
### Script para crear TE fasta

import pandas as pd
from Bio import SeqIO
import sys
import os
from time import gmtime, strftime


# sys.argv[1] == RAL-177
# sys.argv[2] == full-path of genome file of RAL-177

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


# Creo por primera vez el log file

logfile = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2]+ ".log"


### Primero que todo creo un directorio 

dirName = sys.argv[1] + "/" +  sys.argv[2]

if not os.path.exists(dirName):
    os.mkdir(dirName)
else:    
    pass


fo = open(logfile,"w")
log_line = "--- Working with "+ sys.argv[1] + " ---\n"
fo.write(log_line)
log_line = "--- " + strftime("%Y-%m-%d %H:%M:%S", gmtime()) + " ---\n"
fo.write(log_line)
fo.write("--- Creating Genome Folder ---\n")
fo.close()


# Cargo coordenadas de transferencia de TEs y me quedo solo con los genes que estan en cromosomas de interes

TE_annot_file = sys.argv[3]

df_annot_te = pd.read_table(TE_annot_file,sep="\t",names=["seqname","te_start","te_end","teid","te_len","te_strand"])
chromosome_interest = ['4', 'X', '3R', '2L','3L','2R']
df_annot_te = df_annot_te[df_annot_te["seqname"].isin(chromosome_interest)]



genome_file = sys.argv[4]

# Cargo genoma y lo transformo en diccionario

dict_genome = ObtainSeq(genome_file)

# Obtengo la secuencia de los TE

sequence_fasta = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_sequence.fasta"

fo = open(logfile,"a")
log_line = "--- Writing complete TE sequence of "+ sys.argv[2] + " : " + sequence_fasta + " ---\n"
fo.write(log_line)
fo.close()


fo = open(sequence_fasta,"w")

row_list = []
count = 0
TE_n_inners = []
for index,row in df_annot_te.iterrows():
    seq = dict_genome[row.seqname][row.te_start:row.te_end]
    if "N" in seq.upper():
        TE_n_inners.append(row.teid) # Chequeo si hay Ns
    if row.te_strand == "-":
        seq = dict_genome[row.seqname][row.te_start:row.te_end].reverse_complement().seq._data
    else:
        seq = dict_genome[row.seqname][row.te_start:row.te_end].seq._data
    seq = seq.replace("N","") # remplazo las Ns
    fasta = ">" + row.teid + "\n" + seq + "\n"
    fo.write(fasta)
    row_list.append(row.tolist())

fo.close()

# Defino las regiones genomicas a alinear sumando 500 bases a cada extremo del TE

df_annot_te = pd.DataFrame(row_list,columns=df_annot_te.columns.tolist())

flanking_sequence = sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_TE_flanking_regions.fasta"

fo = open(logfile,"a")
log_line = "--- Writing flanking (500 pb) TE sequence of "+ sys.argv[2] + " : " + flanking_sequence + " ---\n"
fo.write(log_line)
fo.close()


fo = open(flanking_sequence,"w")

te_info = []
for index,row in df_annot_te.iterrows():
    pass_filter = "Y"
    contig_boundaries = "N"
    has_N = "N"        
    if row.teid in TE_n_inners:
        has_N = "Y"
    mappeable_size_start = 500
    mappeable_size_end = 500
    if row.te_strand == "-":
        sequence_start = dict_genome[row.seqname][row.te_start-499:row.te_start+1].reverse_complement().seq._data
        sequence_end = dict_genome[row.seqname][row.te_end+1:row.te_end+501].reverse_complement().seq._data
    else:
        sequence_start = dict_genome[row.seqname][row.te_start-501:row.te_start+1].seq._data
        sequence_end = dict_genome[row.seqname][row.te_end+1:row.te_end+499].seq._data
    if "N" in sequence_start.upper():
        if row.te_strand == "+":
            nindex = sequence_start.find("N") # posicion donde esta el indice
            nindex += 100 # Es el tamaño de las N
            mappeable_size_start = len(sequence_start[nindex:]) 
            if mappeable_size_start >= 300:
                sequence_start = sequence_start[nindex:]
            else:
                sequence_start = sequence_start[nindex:]
                contig_boundaries = "Y"
        else:
            nindex = sequence_start.find("N")
            mappeable_size_start = len(sequence_start[:nindex])
            if mappeable_size_start >= 300:
                sequence_start = sequence_start[nindex:]
            else:
                sequence_start = sequence_start[nindex:]
                contig_boundaries = "Y"
    if "N" in sequence_end.upper():
        if row.te_strand == "+":
            nindex = sequence_end.find("N") # posicion donde esta el indice
            mappeable_size_end = len(sequence_end[:nindex]) # Es el tamaño de las N
            if mappeable_size_end >= 300:
                sequence_end = sequence_end[:nindex]
            else:
                sequence_end = sequence_end[:nindex]
                contig_boundaries = "Y"
        else:
            nindex = sequence_end.find("N") + 100  # posicion donde esta el indice
            mappeable_size_end = len(sequence_end[nindex:]) # Es el tamaño de las N
            if mappeable_size_end >= 300:
                sequence_end = sequence_end[:nindex]
            else:
                sequence_end = sequence_end[:nindex]
                contig_boundaries = "Y"
    if mappeable_size_start  == 500 or mappeable_size_end == 500:
        if row.te_strand == "-":
            fasta = ">" + row.teid + "\n" + sequence_end.strip() + sequence_start.strip()  +"\n"
            fo.write(fasta)
        else:
            fasta = ">" + row.teid + "\n" + sequence_start.strip() + sequence_end.strip() +"\n"
            fo.write(fasta)
    else:
        pass_filter = "N"
    single_row = row.tolist() + [pass_filter,has_N,contig_boundaries,mappeable_size_start,mappeable_size_end]
    te_info.append(single_row)
  

fo.close()

TE_sequence_context_info = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_TE_info.tsv"

fo = open(logfile,"a")
log_line = "--- Writing information about the sequence context (inner N, flanking N, etc.) of the TE's "+ sys.argv[1] + " : " + TE_sequence_context_info + " ---\n"
fo.write(log_line)
fo.close()

df_te_info = pd.DataFrame(te_info,columns=["seqname","te_start","te_end","teid","te_len","te_strand","2_transfer","has_n","contig_boundaries","pb_left2map","pb_right2map"])

df_te_info.to_csv(TE_sequence_context_info,sep="\t",header=True,index=False)

fo = open(logfile,"a")
fo.write("--- TE sequence and TE flanking sequence obtained! ---\n")
fo.close()
