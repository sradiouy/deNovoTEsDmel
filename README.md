# deNovoTEsDmel


load the folowwing modules: 
- BEDTools/2.27.1-foss-2018b
-  minimap2/2.11-foss-2016b
-  Python/3.7.2-GCCcore-8.2.0



# Step to follow to transfer TEs and genes:

# Prepare CDS and genes reference fasta

python cds_gene_fasta_preparation.py

#### Set varaible names

strain = base name of a strain (i.e., STO-022)
FastaPath = path to fasta file of genome assemble

### Obtain fasta sequence of TEs from bed file

$ python bedTE2Fasta.py strain FastaPath

### get information of nested and tandem TEs

$ python nested_tandem_TE_classification.py strain

### Gene Transfer 

$ python gene_transfer.py strain FastaPath

### TE Transfer 

$ python TE_transfer.py strain
