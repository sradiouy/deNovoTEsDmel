#!/usr/bin/env python -W ignore::DeprecationWarning
# encoding: utf-8

import pysam
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
import sys
import subprocess 



# Una vez obtenido los CDS y genes de la transferencia a utilizar necesitamos encontrarlos en cada genoma a analizar. Para eso utilizaremos minimap2 en dos modalidades asm5 (para ensamblados genomicos) y splice para los CDS.

##############################################################################################################################
##############################################################################################################################


def runcml(cml):
    """Ejecuta un comando"""
    process = subprocess.Popen(cml.split(),stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=False)
    stdout, stderr = process.communicate()
    stderr = stderr.decode("utf-8")
    print(stderr)
    stdout = stdout.decode("utf-8")
    return stdout


def runminimap2(strain_genome_file,reference_sequence_file,outfile,mode):
    """
    strain_genome_file: es el genoma al cual quiero transferir los genes o cds de la referencia
    reference_sequence_file: son la secuencia de los genes o de los CDS
    mode: genes o cds (archivo fasta del transcripto mas largo)
    """
    if mode == "gene":
        commandLine = "minimap2 -ax asm5 -t 8  -uf -C5 " + strain_genome_file + " " + reference_sequence_file
    elif mode == "cds": # -io ignoro overlap, ignoro -id downstream
        commandLine = "minimap2 -ax splice -t 8 -uf -C5 " + strain_genome_file + " " + reference_sequence_file
    stdout = runcml(commandLine)
    with open(outfile, 'w') as f:
        f.write(stdout)
    return outfile


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


def hamming_distance(s1, s2):
    """
    Determino las diferencias entre dos secuencias
    """
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def create_synteny_refrence(dataframe):
    synteny = []
    for index,row in dataframe.iterrows():
        if index == 0:
            prev_gene = "-"
            prev_chr = "-"
            post_gene = dataframe.loc[index+1].gene
            post_chr = dataframe.loc[index+1].seqname
        elif index == len(dataframe) -1:
            post_gene = "-"
            post_chr = "-"
            prev_gene = dataframe.loc[index-1].gene
            prev_chr = dataframe.loc[index-1].seqname
        else:
            post_gene = dataframe.loc[index+1].gene
            prev_gene = dataframe.loc[index-1].gene
            post_chr = dataframe.loc[index+1].seqname
            prev_chr = dataframe.loc[index-1].seqname
        if prev_chr != row.seqname:
            prev_gene = "-"
        if post_chr != row.seqname:
            post_gene = "-"
        synteny.append(row.tolist() + [prev_gene,post_gene])
    header = dataframe.columns.tolist() + ["prev_gene","post_gene"]
    dataframe = pd.DataFrame(synteny,columns=header)
    return dataframe


def create_synteny(dataframe):
    synteny = []
    for index,row in dataframe.iterrows():
        if index == 0:
            prev_gene = "-"
            prev_chr = "-"
            post_gene = dataframe.loc[index+1].gene
            post_chr = dataframe.loc[index+1].chromosome
        elif index == len(dataframe) -1:
            post_gene = "-"
            post_chr = "-"
            prev_gene = dataframe.loc[index-1].gene
            prev_chr = dataframe.loc[index-1].chromosome
        else:
            post_gene = dataframe.loc[index+1].gene
            prev_gene = dataframe.loc[index-1].gene
            post_chr = dataframe.loc[index+1].chromosome
            prev_chr = dataframe.loc[index-1].chromosome
        if prev_chr != row.chromosome:
            prev_gene = "-"
        if post_chr != row.chromosome:
            post_gene = "-"
        synteny.append(row.tolist() + [prev_gene,post_gene])
    header = dataframe.columns.tolist() + ["prev_gene","post_gene"]
    dataframe = pd.DataFrame(synteny,columns=header)
    return dataframe

def compare_strain_refrence_synteny(synteny_strain,synteny_reference):
    list_rows = []
    for index,row in synteny_strain.iterrows():
        prev_gene,post_gene = synteny_reference[synteny_reference["gene"] == row.gene][["prev_gene","post_gene"]].values[0].tolist()
        if prev_gene == row.prev_gene:
            synteny_prev = "Y" # Si para el gen en cuestion el gen previo es el mismo que en la transferencia ponemos Y
        else:
            synteny_prev = "N" # Si para el gen en cuestion el gen previo es distinto que en la transferencia ponemos N
        if post_gene == row.post_gene:
            synteny_post = "Y" # Si para el gen en cuestion el gen posterior es el mismo que en la transferencia ponemos Y
        else:
            synteny_post = "N" # Si para el gen en cuestion el gen posterior es el mismo que en la transferencia ponemos N
        diff_ref = abs(synteny_reference[synteny_reference["gene"] == row.gene].transfer_length.values[0] - row.length) # Calculamos la diferencia de tamaño entre el CDS del strain y de la referencia
        if synteny_reference[synteny_reference["gene"] == row.gene].transfer_length.values[0] > row.length:
            bigger = "Y" # Si el strain es mas grande
        elif synteny_reference[synteny_reference["gene"] == row.gene].transfer_length.values[0] == row.length:
            bigger = "-" # Si el strain es igual a la referencia
        else:
            bigger = "N" # Si el strain es mas chico
        transfer_start = synteny_reference[synteny_reference["gene"] == row.gene].transfer_start.values[0]
        diff_start_pos = abs(transfer_start-row.start)
        list_rows.append(row.tolist() + [synteny_prev,synteny_post,diff_ref,bigger,diff_start_pos])
    # Guardo la info en una nueva tabla
    header = synteny_strain.columns.tolist() + ["synteny_prev","synteny_post","diff_ref","bigger","diff_start_pos"]
    synteny_strain = pd.DataFrame(list_rows,columns=header)
    return synteny_strain

### INICIO

### Abro el log 

logfile = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2]+ ".log"


fo = open(logfile,"a")
fo.write("--- Gene transfer started ---\n")
fo.close()


# Cargo secuencia genomica en forma de un diccionario
genome = ObtainSeq(sys.argv[3])


# Cargo coordenadas de transferencia de CDS y me quedo solo con los genes que estan en cromosomas de interes

transfer_file = "Reference_primary_cds_transfer.bed"

df_transfer = pd.read_table(transfer_file,sep="\t",names=["seqname","transfer_start","transfer_end","gene","transfer_frame","transfer_strand"])

df_transfer["transfer_length"] = abs(df_transfer["transfer_end"] - df_transfer["transfer_start"])

chromosome_interest = ['4', 'X', '3R', '2L','3L','2R']

df_transfer = df_transfer[df_transfer["seqname"].isin(chromosome_interest)]

# Me quedo unicamente con los genes de interes (que son los que estan en el )

genes_interest = df_transfer.gene.tolist()


#### Analisis con CDS

cds_sequence_file = "Reference_primary_cds.fasta"

cdssamfile =  sys.argv[1] + "/" + sys.argv[2] + "/cds2" + sys.argv[2] + ".sam"


fo = open(logfile,"a")
log_line = "--- Number of gene considered for transfer: " + str(len(genes_interest)) + "  ---\n"
fo.write(log_line)
fo.close()


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for CDS sequences ---\n")
fo.close()

runminimap2(sys.argv[3],cds_sequence_file,cdssamfile,"cds")


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for CDS sequences: Done ---\n")
fo.close()

## Defino todos los reads existentes y los llevo a un dataframe

samfilete = pysam.AlignmentFile(cdssamfile, "r")

list_te = []
for read in samfilete.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    list_te.append([read.query_name,read.reference_start,read.reference_end,read.is_unmapped])

samfilete.close()

df_to_define_type_of_mapping = pd.DataFrame(list_te,columns=["gene","start","end","unmapped"])


rows = []
for name,grp in df_to_define_type_of_mapping.groupby("gene"):
    if grp.start.values[0] == -1:
        rows.append([sys.argv[1],name,0,"NA"])
    else:
        if len(grp) == 1:
            rows.append([sys.argv[1],name,1,"unique"])
        else:
            rows.append([sys.argv[1],name,len(grp),"MM"])
    
df_to_define_type_of_mapping = pd.DataFrame(rows,columns=["genome","gene","mmcount","type"])


mapping_information_file =  sys.argv[1] + "/"  + sys.argv[2] + "/" + sys.argv[2] + "_CDS_mapping_information.tsv"


df_to_define_type_of_mapping.to_csv(mapping_information_file,sep="\t",index=False,header=True)


fo = open(logfile,"a")
log_line = "--- Mapping information of CDS obtained: " + mapping_information_file+" ---\n"
fo.write(log_line)
fo.close()


# Cargo el resultado de minimap2 (un archivo SAM) a un objeto pysam, para poder iterar cada read del resultado


fo = open(logfile,"a")
fo.write("--- Transfer of CDS position is starting ---\n")
fo.close()

samfilecds = pysam.AlignmentFile(cdssamfile, "r")

# Itero el resultado

rows = [] # genero una lista vacia para ir guardando la informacion de cada mapeo
duplicated_counts = {} # sera un contador para establecer cuantos mapeos primarios me caen en el mismo lugar
missedcds = [] # para saber cuantos cds no son encontrados
shortmissedcds = [] # para saber cuantos de los perdidos son chicos y ver si lo podemos recuperar
all_genes = [] # para contar el numero inicial de genes analizados

# establezco el header de la tabla a crear
header = ["chromosome","start","end","gene","frame","strand","NM","nm_start","nm_end","length","mapq","cigar","flag","read_id"]

for read in samfilecds.fetch(): # obtengo la informacion por cada read. Convierto el read en un clase con sus metodos
    gene = read.query_name # me quedo con el nombre del read (en este caso es el gen del cual se origino el CDS)
    if gene not in genes_interest: # Si el gen no esta entre los genes de interes lo ignoro
        pass
    else: 
        NM = 0 # inicializo el numero de missmatch del read
        nm_start = 0 # inicializo el numero de missmatch del read al inicio de la secuencia (20 pb)
        nm_end = 0 # inicializo el numero de missmatch del read al inicio de la secuencia (20 pb)
        chromosome = read.reference_name # establezco el nombre del cromosoma
        start = read.reference_start # establezco la posicion de inicio en el genoma de referencia (coordenada de mapeo)
        end = read.reference_end # establezco la posicion de fin en el genoma de referencia (coordenada de mapeo + CIGAR)
        all_genes.append(gene) # guardo el nombre del gen
        if read.is_reverse: # establezco el strand
            strand = "-"
        else:
            strand = "+"
        length = read.reference_length # Me quedo con el largo de la referencia
        mappingquality = read.mapping_quality # Me quedo con la calidad de mapeo
        cigar = read.cigarstring # guardo el CIGAR
        if start == -1 or read.is_unmapped: # sino mapeo el start se setea como -1
            missedcds.append(gene) # por lo tanto lo guardo en missedcds
            if read.query_length <= 150: # Si es pequeño puede ser que el algoritmo falle
                shortmissedcds.append(gene) # por lo tanto lo guardo en shortmissedcds
        elif read.flag == 272 or read.flag == 256: # Si no es alineamiento primario o suplementario, espero a conseguirlo
            pass
        else:
            read_id = chromosome + "_" + str(start) +  "_" + str(end) # genero coordenado de mapeo unica (para aquellos genes MMappping)
            if gene in duplicated_counts.keys(): 
                duplicated_counts[read_id] += 1 # Si existe sumo 1
            else:
                duplicated_counts[read_id] = 1 # Si no existe lo inicio 1
            nm_start = hamming_distance(genome[chromosome][start:end].seq._data[:20], read.query_alignment_sequence[:20]) # Comparo las primeras 20 bases de la secuencia de referencia con las primeras de la referencia y determino el n de missmatch
            nm_end = hamming_distance(genome[chromosome][start:end].seq._data[-20:], read.query_alignment_sequence[-20:])
            # Comparo las ultimas 20 bases de la secuencia de referencia con las primeras de la referencia y determino el n de missmatch
            NM = [x[1] for x in read.tags if "NM" in x] # guardo el numero de missmatch de el alineamiento entero
            rows.append([chromosome,start,end,gene,".",strand,NM[0],nm_start,nm_end,length,mappingquality,cigar,read.flag,read_id]) # guardo la info en una tabla


samfilecds.close()

## Genes perdidos 


fo = open(logfile,"a")
log_line = "--- Missed CDS sequences: " + str(len(missedcds)) + " ---\n"
fo.write(log_line)
fo.close()

## genes perdidos y pequeños


fo = open(logfile,"a")
log_line = "--- Missed CDS sequences of short size (< 150): " + str(len(shortmissedcds)) + " ---\n"
fo.write(log_line)
fo.close()

all_genes = list(set(all_genes)) # me quedo con genes unicos

dfcds = pd.DataFrame(rows,columns=header)

fo = open(logfile,"a")
log_line = "--- CDS mapped: " + str(len(dfcds.gene.unique())) + " ---\n"
fo.write(log_line)
fo.close()

# Cuento cuantos multimmapping
mmappinggenes = []
problem_regions = 0
for name,grp in dfcds.groupby("read_id"):
    if len(grp) >= 2:
        problem_regions += 1
        mmappinggenes += grp.gene.tolist()


fo = open(logfile,"a")
log_line = "--- Multimapping regions: " + str(problem_regions) + " ---\n"
fo.write(log_line)
fo.close()


dfcds = dfcds.drop(['read_id'], axis=1) # elimino la columna de read_id

### Trato de recuperar los genes perdidos que son de pequeño tamaño

fo = open(logfile,"a")
fo.write("--- Recovering CDS sequence of short size ---\n")
fo.close()

# Para eso la estrategia es ver si mapeando los genes con minimap2 somos capaces de determinar su ubicacion. Luego transformaremos las coordenadas de inicio y fin del mapeo a las CDS localizando el start y stop codon dentro del gen.

# Para eso utilizaremos Biopython para crear un diccionario de fastas de cada CDS

cdsfasta = "Reference_primary_cds.fasta"

cdsfasta = open(cdsfasta)
cdsfastadict = SeqIO.to_dict(SeqIO.parse(cdsfasta, "fasta"))


# Cargo el resultado de minimap2 de genes (un archivo SAM) a un objeto pysam, para poder iterar cada read del resultado


gene_sequence_file = "/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/Reference_genes.fasta"

genesamfile =  sys.argv[1] + "/" + sys.argv[2] + "/gene2" + sys.argv[2] + ".sam"

fo = open(logfile,"a")
fo.write("--- Running Minimap2 for gene sequences ---\n")
fo.close()

runminimap2(sys.argv[3],gene_sequence_file,genesamfile,"gene")


fo = open(logfile,"a")
fo.write("--- Running Minimap2 for gene sequences: Done ---\n")
fo.close()

samfilegenes = pysam.AlignmentFile(genesamfile, "r")

rows = [] # genero una lista vacia para ir guardando la informacion de cada mapeo
missedgenes = [] # para saber cuantos genes no son encontrados

# establezco el header de la tabla a crear

header = ["chromosome","start","end","gene","frame","strand","NM","nm_start","nm_end","length","mapq","cigar","flag"]

# Las variables son las mismas a lo hecho anteriormente

for read in samfilegenes.fetch():
    NM = 0
    nm_start = 0
    nm_end = 0
    chromosome = read.reference_name
    start = read.reference_start
    end = read.reference_end
    gene = read.query_name
    if read.is_reverse:
        strand = "-"
    else:
        strand = "+"
    length = read.reference_length
    mappingquality = read.mapping_quality
    cigar = read.cigarstring
    if start == -1:
        missedgenes.append(gene)
    elif read.flag == 272 or read.flag == 256 or read.flag == 2048 or read.flag == 2064: # Elimino tambien los alineamientos suplementarios
        pass
    else:
        if gene in shortmissedcds: # solo analizo los genes (reads) que estan entre los perdidos
            ref_start = 0
            ref_end = 0
            nm_start = hamming_distance(genome[chromosome][start:end].seq._data[:20], read.query_alignment_sequence[:20])
            nm_end = hamming_distance(genome[chromosome][start:end].seq._data[-20:], read.query_alignment_sequence[-20:])
            # Determino el inicio y fin del CDS en el gen, con un match perfecto de 30 nucleotidos para estar completamente seguros, si no es asi lo paso.
            if strand == "-":
                sequence = read.query_alignment_sequence
                sequence_start = cdsfastadict[gene].seq.reverse_complement()._data[:30] # determino el inicio y hago la reversa complementaria ya que mapeo en la hebra -
                if sequence.count(sequence_start) == 1: # cuento si esa secuencia esta presente el la region mapeada
                    ref_start = start + sequence.find(sequence_start) # si esta busco cual es la posicion del indice (donde se encuentra) y genero un nuevo start
                else:
                    ref_start = 0 # si no esta no marco ninguna coordenada nueva
                sequence_end = cdsfastadict[gene].seq.reverse_complement()._data[-30:] # determino el fin y hago la reversa complementaria ya que mapeo en la hebra -. Luego aplico la misma estrategia.
                if sequence.count(sequence_end) == 1:
                    ref_end = start + (sequence.find(sequence_end) + 30)
            else: # En el caso que la hebra sea + la estrategia es la misma sin la necesidad de hacer reveresa complementaria
                sequence = read.query_alignment_sequence
                sequence_start = cdsfastadict[gene].seq._data[:30]
                if sequence.count(sequence_start) == 1:
                    ref_start = start + sequence.find(sequence_start)
                else:
                    ref_start = 0
                sequence_end = cdsfastadict[gene].seq._data[-30:]
                if sequence.count(sequence_end) == 1:
                    ref_end = start + (sequence.find(sequence_end) + 30)    
                else:
                    ref_end = 0
            if ref_start == 0 or ref_end == 0:
                pass
            else:
                start = ref_start
                end = ref_end
                NM = [x[1] for x in read.tags if "NM" in x][0]
                rows.append([chromosome,start,end,gene,".",strand,NM,nm_start,nm_end,length,mappingquality,cigar,read.flag])

samfilegenes.close()

dfgenes = pd.DataFrame(rows,columns=header)

## Genes recuperados 

fo = open(logfile,"a")
log_line = "--- " + str(len(dfgenes)) + " CDS recovered ---\n"
fo.write(log_line)
fo.close()

# Agrego los CDS recuperados al dataframe de CDS

dfcds = pd.concat([dfcds,dfgenes])

fo = open(logfile,"a")
log_line = "--- In total there were " + str(len(dfcds.gene.unique())) + " genes/CDS transferred ---\n"
fo.write(log_line)
fo.close()

# Una vez obtenido los resultados de los mapeos, queremos ver como quedo la sintenia, definida unicamente por los genes que estan antes y despues del gen en estudio. Es solo mirando coordenadas sin importar la hebra. 


# Ordeno por cromosoma posicion de incio y fin

dfcds = dfcds.sort_values(["chromosome","start","end"])
dfcds = dfcds.reset_index(drop=True)

# Establezco la sintenia

fo = open(logfile,"a")
fo.write("--- Gene synteny comparision between strain-reference ---\n")
fo.close()

dfcds = create_synteny(dfcds)

# Como la idea es comporbar como se comporta la sintenia necesito determinar la sintenia del archivo de CDS de transferencia

df_transfer_original = df_transfer # creo backup

# Me quedo solamente con los genes que pude transferir
df_transfer = df_transfer[df_transfer["gene"].isin(list(dfcds.gene.unique()))]

# Y ordeno el dataframe por cromosoma y coordenadas

df_transfer = df_transfer.sort_values(["seqname","transfer_start","transfer_end"])
df_transfer = df_transfer.reset_index(drop=True)

## De la misma forma que antes establezco la sintenia y creo la nueva tabla con esa info.

df_transfer = create_synteny_refrence(df_transfer)

# Una vez establecida ambas sintenias compararemos ambos casos para ver como se comportan tanto el gen previo como el post entre la referencia y el strain (veo si se mantiene)

df_cds_synteny = compare_strain_refrence_synteny(dfcds,df_transfer)

# Como en el mapeo existen mapeos suplementarios el objetivo es determinar cual es el real y para eso veremos cual de ellos mantiene la sintenia y cual esta mas cerca a la coordenada de inicio de la referencia

df_cds_synteny = df_cds_synteny.sort_values(["gene","synteny_prev","synteny_post","diff_start_pos"],ascending = [False,False,True,True])

# eliminamos los duplicados 

df_cds_synteny = df_cds_synteny.drop_duplicates("gene",keep="first")

removed_mm = len(dfcds) - len(df_cds_synteny)


fo = open(logfile,"a")
log_line = "--- In this step we eliminate duplicated genes based in the synteny, we keep the position were the synteny is mantained, so removing probably ture duplications. In total we removed  " + str(removed_mm) + " genes/CDS duplications ---\n"
fo.write(log_line)
fo.close()


# Una vez eliminado estos duplicados procederemos a calcular nuevamente la sintenia, por lo que eliminamos la info de gen previo y posterior


fo = open(logfile,"a")
fo.write("--- Redone of the gene synteny comparision between strain-reference with the new set of CDS ---\n")
fo.close()


df_cds_synteny = df_cds_synteny[['chromosome', 'start', 'end', 'gene', 'frame', 'strand', 'NM','nm_start', 'nm_end', 'length', 'mapq', 'cigar', 'flag']]

# Ordenamos nuevamente y reseteamos el indice (ordenar siempre cambia el orden del indice)

df_cds_synteny = df_cds_synteny.sort_values(["chromosome","start","end"])
df_cds_synteny = df_cds_synteny.reset_index(drop=True)

df_cds_synteny = create_synteny(df_cds_synteny)

# Una vez establecida ambas sintenias volvemos a comparar ambos casos para ver como se comportan tanto el gen previo como el post entre la referencia y el strain (veo si se mantiene)

df_cds_synteny = compare_strain_refrence_synteny(df_cds_synteny,df_transfer)

# Calculo el tamaño de la CDS en la referencia y lo agrego a la tabla 

df_cds_synteny["reference_size"] = [df_transfer[df_transfer["gene"] == gene].transfer_length.values[0] for gene in df_cds_synteny.gene.tolist()]

# Determino en porcentaje cuan diferente es el tamaño entre el CDS del strain y de la referencia 

df_cds_synteny["percent_size"] = (df_cds_synteny["reference_size"] / df_cds_synteny["length"])*100
df_cds_synteny["percent_size"] = df_cds_synteny["percent_size"].astype(int)

## Veo cuantos genes en total se perdieron

len_missed_genes = len(all_genes) - len(df_cds_synteny)


fo = open(logfile,"a")
log_line = "--- Total missed genes: " + str(len_missed_genes)+ " ---\n"
fo.write(log_line)
fo.close()


missed_genes = [ gene for gene in all_genes if gene not in df_cds_synteny.gene.tolist()]

# Determino cuantos gene problematicos hay. Un gen problematico es aquel que el tamaño es mayor o menor al 105 o 95 % respectivamente y ademas no mantuvo la sintenia en alguno de los dos genes

len_problematic_genes = len(df_cds_synteny[((df_cds_synteny["percent_size"] > 105) | (df_cds_synteny["percent_size"] < 95))& ((df_cds_synteny.synteny_prev == "N") | (df_cds_synteny.synteny_post == 'N'))])


fo = open(logfile,"a")
log_line = "--- Remove of problematics genes (percent size > 105 or < 95 and not loss of synteny): "+ str(len_problematic_genes)+" genes removed from the bed transfer file ---\n"
fo.write(log_line)
fo.close()


problematic_genes = df_cds_synteny[((df_cds_synteny["percent_size"] > 105) | (df_cds_synteny["percent_size"] < 95))& ((df_cds_synteny.synteny_prev == "N") | (df_cds_synteny.synteny_post == 'N'))].gene.tolist()


missed_genes_file = "/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Results/" + sys.argv[1] + "/" + sys.argv[1] +"_gene_absent.tsv"

fo = open(logfile,"a")
fo.write("--- Writing missed/problematic gene file ---\n")
fo.close()

# Creo una lista de los genes a no incluir y la razon de ello.
absent = missed_genes +  problematic_genes
absent_info = ["Not mapped"]*len_missed_genes + ["No perfect synteny or length"]*len_problematic_genes 

# genero un dataframe y guardo la info

df_missing_genes = pd.DataFrame(
    {'gene': absent,
     'reason': absent_info})

df_missing_genes.to_csv(missed_genes_file,sep="\t",index=False,header=True)


### Chequeo problematicos, inversiones, MM y problemas relacionados con cromsomas


## Marco los genes problematicos por tamaño 
df_cds_synteny["Problematic"] = ["Y" if gene in problematic_genes else "N" for gene in df_cds_synteny.gene.tolist()]


### Marco los MM genes y aquellos vecinos a ellos


fo = open(logfile,"a")
fo.write("--- Adding multimapping information ---\n")
fo.close()


new_list = []
for index,row in df_cds_synteny.iterrows():
    mmcount = df_to_define_type_of_mapping[df_to_define_type_of_mapping["gene"]==row.gene].mmcount.values[0]
    try:
        mmcount_prev = df_to_define_type_of_mapping[df_to_define_type_of_mapping["gene"]==row.prev_gene].mmcount.values[0]
    except:
        mmcount_prev = 0
    try:
        mmcount_pos = df_to_define_type_of_mapping[df_to_define_type_of_mapping["gene"]==row.post_gene].mmcount.values[0]
    except:
        mmcount_pos = 0 
    if mmcount > 1: 
        mmapping = "MM" # gen Multimapping
    elif mmcount_pos > 1 or mmcount_prev > 1:
        mmapping = "RMM" # relacionado a un gen Multimapping
    else:
        mmapping = "NMM" # No esta relacionado a Multimapping
    new_list.append(row.tolist() + [mmapping])

header = df_cds_synteny.columns.tolist() + ["mmapping"]

df_cds_synteny_mm = pd.DataFrame(new_list,columns=header)


#### MARCO LAS INVERSIONES 

fo = open(logfile,"a")
fo.write("--- Adding inversion information ---\n")
fo.close()

inverted_genes = []
for index,row in df_cds_synteny_mm.iterrows():
    if row.strand != df_transfer[df_transfer["gene"] == row.gene].transfer_strand.values[0]:
        if row.gene not in inverted_genes:
            inverted_genes.append(row.gene)

# Determino y marco las inversiones y los genes posiblemente afectados por ella

new_list = []
for index,row in df_cds_synteny_mm.iterrows():
    if row.gene in inverted_genes: 
        inversion = "I" # gen invertido
    elif row.post_gene in inverted_genes or row.prev_gene in inverted_genes:
        inversion = "RI" # relacionado a un gen invertido
    else:
        inversion = "NI" # No esta relacionado
    new_list.append(row.tolist() + [inversion])

header = df_cds_synteny_mm.columns.tolist() + ["inversion"]

df_cds_synteny_mm_inversion = pd.DataFrame(new_list,columns=header)


#### Marco los genes que tiene diferencia en el cromosoma entre la referencia y el strain

fo = open(logfile,"a")
fo.write("--- Adding chromosome translocation information ---\n")
fo.close()

new_list = []
for index,row in df_cds_synteny_mm_inversion.iterrows():
    if row.chromosome != df_transfer[df_transfer["gene"]== row.gene].seqname.values[0]:
        chromosome_problem = "C" # diferente coromosoma que referencia
    else:
        chromosome_problem = "NC" # No esta relacionado
    new_list.append(row.tolist() + [chromosome_problem])
 

header = df_cds_synteny_mm_inversion.columns.tolist() + ["chromosome_problem"]

df_cds_synteny_mm_inversion_c = pd.DataFrame(new_list,columns=header)

## Veo los genes donde no coincide la sintenia de forma completa

# Marco problemas de sintenia

fo = open(logfile,"a")
fo.write("--- Adding synteny information ---\n")
fo.close()

df_cds_synteny_mm_inversion_c_s = df_cds_synteny_mm_inversion_c

synteny_check = df_cds_synteny_mm_inversion_c['synteny_prev'] + df_cds_synteny_mm_inversion_c['synteny_post']

# Si los genes previos y posteriores no coinciden con la referencia pongo como si tuviera un problema

df_cds_synteny_mm_inversion_c_s["synteny_problem"] = ["N" if x == "YY" else "Y" for x in synteny_check.tolist()]

# Guardo la tabla

gene_information_file = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer_info.tsv"

fo = open(logfile,"a")
log_line = "--- Writing gene information table: " + gene_information_file + " ---\n"
fo.write(log_line)
fo.close()

df_cds_synteny_mm_inversion_c_s.to_csv(gene_information_file,sep="\t",index=False,header=True)


# Creo un bed con los genes sin los genes descartados y marcando problema de sintenia, mmapping, coromosoma e inversiones 

gene_bed_file = sys.argv[1] + "/" + sys.argv[2] + "/" + sys.argv[2] + "_gene_transfer.bed"


fo = open(logfile,"a")
log_line = "--- Writing BED file of genes transferred to the strain genome (removing missed and problematics genes): " + gene_bed_file + " ---\n"
fo.write(log_line)
fo.close()

present_gene = df_cds_synteny_mm_inversion_c_s[~df_cds_synteny_mm_inversion_c.gene.isin(df_missing_genes.gene.tolist())]

present_gene = present_gene[['chromosome', 'start', 'end', 'gene', 'frame', 'strand',"mmapping","inversion","chromosome_problem","synteny_problem"]]

present_gene.to_csv(gene_bed_file,sep="\t",index=False,header=False)

fo = open(logfile,"a")
log_line = "--- In total there were " + str(len(present_gene)) + " gene transferred ---\n"
fo.write(log_line)
fo.close()

# Creo un bed de la referencia utilizando unicamente los genes que tranfiero:

df_transfer_subset = df_transfer[df_transfer["gene"].isin(present_gene.gene.tolist())]

df_transfer_subset = df_transfer_subset.sort_values(["seqname","transfer_start","transfer_end"])
df_transfer_subset = df_transfer_subset.reset_index(drop=True)

df_transfer_subset = df_transfer_subset[['seqname', 'transfer_start', 'transfer_end', 'gene', 'transfer_frame','transfer_strand']]

reference_gene_bed_file = "/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Results/"+sys.argv[1]+"/Reference_gene_for_"+sys.argv[1]+"_transfer.bed"

df_transfer_subset.to_csv(reference_gene_bed_file,sep="\t",header=False,index=False)

fo = open(logfile,"a")
log_line = "--- Writing BED file assocciated with the reference genome containing only the genes annotated in the strain genome (removing missed and problematics genes): " + reference_gene_bed_file + " ---\n"
fo.write(log_line)
fo.close()

fo = open(logfile,"a")
fo.write("--- Gene transfer done ---\n")
fo.close()
