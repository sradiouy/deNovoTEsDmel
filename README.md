# TE TRANSFER


#### Programs versions used:

-  BEDTools/2.27.1-foss-2018b
-  minimap2/2.11-foss-2016b
-  Python/3.7.2-GCCcore-8.2.0


# Step to follow to transfer TEs and genes:

#### Set varaible names

- wdir = working directory
- strain = base name of a strain (i.e., STO-022)
- FastaPath = full path to fasta file of genome assemble
- TEAnnotation = full path to TE annotation file in bed format

### Prepare CDS and genes reference fasta

#### Download the following files from flyBase and change accordingly in the following script:

- dmel-all-r(VERSION).gtf'
- dmel-all-gene-r(VERSION).fasta
- dmel-all-CDS-r(VERSION).fasta

python cds_gene_fasta_preparation.py wdir

### Obtain fasta sequence of TEs from bed file

$ python bedTE2Fasta.py wdir strain TEAnnotation FastaPath

### get information of nested and tandem TEs

$ python nested_tandem_TE_classification.py strain

### Gene Transfer 

$ python gene_transfer.py strain FastaPath

### TE Transfer 

$ python TE_transfer.py strain
