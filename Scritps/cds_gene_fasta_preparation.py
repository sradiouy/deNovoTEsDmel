import pysam
from collections import defaultdict
import pandas as pd
from Bio import SeqIO


# El primer paso de la estrategia consiste en la preparacion de los archivos fasta tanto de CDS como
# de genes.  

# Por lo tanto cargaremos el archivo fasta al sistema. Y lo usaremos inicialmente para detectar cuales CDS 
# cominezan con un ATG y cuales terminan con un STOP

cds_fasta = '/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/dmel-all-CDS-r6.31.fasta'

transcripts = [] # Defino una lista vacia donde guardaremos cada CDS (Transcripto y CDS lo uso de forma equivalente)
stopcodons = ["TAA","TAG","TGA"] # Defino los codones STOP
with open(cds_fasta,"r") as fi: # Comienzo a leer linea a linea el archivo
    prev_line = ""
    for line in fi:
        if ">" in line: # Cuando es un header entro aqui
            if prev_line != "": # Como quiero que este atrasado una linea hago que la primera vez no entre.
                if prev_line.strip()[-3:].upper() in stopcodons: # Luego determino si el codon es o no stop
                    stop = "Y"
                else:
                    stop = "N"
                transcripts.append([name,chromosome,gene,transcript,length,atg,stop]) # agrgo informacion a la lista de transcriptos nascientes
            count = 0
            gene,transcript = line.strip("").split(";")[-4].strip().split(",") # parseo el header para quedarme con el nombre del gen y del transcripto
            gene = gene.split("parent=")[1] # parseo aun mas el nombre del gen
            length = int(line.strip("").split(";")[-5].strip().split("length=")[1]) # me quedo con el largo
            name = line.split(" ")[0].replace(">","") # guardo el nombre
            chromosome = line.strip("").split(";")[1].strip().split(":")[0].replace("loc=","") # parseo el header
        else: # Si es secuencia
            prev_line = line # guardo la linea para luego verificar el codon stop
            count += 1
            if count == 1:
                if line[:3].upper() == "ATG": # chequeo si comienza con ATG
                    atg = "Y"
                else:
                    atg = "N"
                    

             
dfcds = pd.DataFrame(transcripts,columns=["name","chromosome","gene","transcripts","length","atg","stop"])  # guardo esa info en un dataframe      

# Me quedo solo con los brazos mayores

dfcds = dfcds[dfcds["chromosome"].isin(['2L', '2R', '3L', '3R', 'X'])]

# Luego la idea es quedanos con un solo transcripto por gen. Para eso agrupo por gen y la idea es quedarnos con el transcripto cuyo CDS sea mas largo. Y si el largo es igual me quedo con los que tenga inicio y stop. 

list_rows = []
for name,grp in dfcds.groupby("gene"):
    subset_grp = grp.sort_values(["length","atg","stop","transcripts"],ascending=[False,False,False,True])
    count = 0
    for index,row in grp.iterrows():
        count += 1
        if count == 1:
            rows = row.tolist() + ["Y"]*2
        else:
            rows = row.tolist() + ["Y","N"]
        list_rows.append(rows)


dfcds_info = pd.DataFrame(list_rows,columns=["name","chromosome","gene","transcript","length","atg","stop","gene_included","transcript_included"])         

dfcds_info.to_csv("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/Reference_primary_cds_info.tsv",sep="\t",index=False,header=True) 


#### Finalmente escribo el nuevo fasta con solamente con los CDS derivados del transcripto de interes.

with open(cds_fasta,"r") as fi, open("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/Reference_primary_cds.fasta","w") as fo:
    for line in fi:
        if ">" in line:
            flag = False
            gene,transcript = line.strip("").split(";")[-4].strip().split(",")
            gene = gene.split("parent=")[1]
            if dfcds_info[dfcds_info["transcript"]==transcript]["transcript_included"].values == "Y":
                fo.write(">" + gene + "\n")
                flag = True
        elif flag:
            fo.write(line)
            
 
# Aplico la misma estrategia con genes, quedandome únicamente con aquellos genes que tengan un representante en el fasta de CDS.  
            
genes = list(dfcds_info.gene.unique()) # Determino cuales son esos genes

genes_fasta = '/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/dmel-all-gene-r6.31.fasta' # establezo el archivo 

# Creo el fasta de los genes seleccionados. 

with open(genes_fasta,"r") as fi, open("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/genes.fasta","w") as fo:
    for line in fi:
        if ">" in line:
            gene_flag = False
            gene = line.split(" ")[0].replace(">","")
            if gene in genes:
                gene_flag = True
                fo.write(">" + gene + "\n")
        elif gene_flag:
            fo.write(line)
 
 
### Luego genero un bed madre para los CDS y genes seleccionados a partir del gtf
            
gff_file = '/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/dmel-all-r6.31.gtf'

header = ["seqname","source","feature","start","end","score","strand","frame","attribute"]

df_gff = pd.read_table(gff_file,sep="\t",header=None,names=header,comment="#")


### Genero bed de genes usando como feature el gen 

df_genes = df_gff[df_gff["feature"] == "gene"]


genes_rows = []
for index,row in df_genes.iterrows():
    gene = row.attribute.split(";")[0].split(' ')[1].replace('"',"")
    if gene in genes:
        single = [row.seqname,row.start,row.end,gene,".",row.strand]
        genes_rows.append(single)

pd.DataFrame(genes_rows).to_csv("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/Reference_genes_transfer.bed",sep="\t",header=False,index=False)

### Genero bed de CDSs usando el codon start y codon stop como marcadores del inicio y fin del transcripto

## Primero me quedo con los STOP y START de los CDS de interes

df_cds = df_gff[(df_gff["feature"] == "start_codon") | (df_gff["feature"] == "stop_codon")]
df_cds["transcripts"] = df_cds.attribute.str.split(";").str[2].str.split(" ").str[-1].str.replace('"',"")
df_cds["gene"] =df_cds.attribute.str.split(";").str[0].str.split(" ").str[-1].str.replace('"',"")
df_cds = df_cds[["seqname","feature","start","end","strand","transcripts","gene"]]

cdss = list(dfcds_info[dfcds_info["transcript_included"] == "Y"].transcript.unique())
# Guardo los nombres de ellos


# Establezco las cordenadas de inicio y fin de cada CDS y lo transformo en un bed

cdss_rows = []

for singlecds in cdss:
    subdf = df_cds[df_cds["transcripts"] == singlecds]
    for index, row in subdf.iterrows():
        gene = row.gene
        if row.feature == "start_codon":
            if row.strand == "+":
                start = row.start
            else:
                end = row.start
        if row.feature == "stop_codon":
            if row.strand == "+":
                end = row.end
            else:
                start = row.end
    cdss_rows.append([row.seqname,start,end,gene,".",row.strand])
    

header = ["seqname","start","end","gene","frame","strand"]

df_cds_bed = pd.DataFrame(cdss_rows,columns=header)


chromosome_interest = ['4', 'X', '3R', '2L','3L','2R']

df_cds_bed = df_cds_bed[df_cds_bed["seqname"].isin(chromosome_interest)]

df_cds_bed = df_cds_bed.sort_values(["seqname","start","end"])
df_cds_bed = df_cds_bed.reset_index(drop=True)
   
df_cds_bed.to_csv("/homes/users/sradio/scratch/eQTL_Dros/TEs_genomes_annotation/minimap2/Annotation_Files/Reference_primary_cds_transfer.bed",sep="\t",header=False,index=False)
