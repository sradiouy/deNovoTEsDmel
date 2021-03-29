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
- RefrencePath =  full path to Reference genome in fasta format
- genes2strain = sam file transfer of genes from reference genomes to strains
- cds2strain = sam file transfer of cds from reference genomes to strains

### Prepare CDS and genes reference fasta

#### Download the following files from flyBase and change accordingly in the following script:

- dmel-all-r(VERSION).gtf'
- dmel-all-gene-r(VERSION).fasta
- dmel-all-CDS-r(VERSION).fasta

python cds_gene_fasta_preparation.py wdir

### TE fasta sequence

$ python bedTE2Fasta.py wdir strain TEAnnotation FastaPath

### Nested and tandem TEs

$ python nested_tandem_TE_classification.py wdir strain TEAnnotation

### Gene transfer

#### Run Minimap2 in splice-aware mode, and without it (Reference_genes.fasta and Reference_primary_cds.fasta: produced by cds_gene_fasta_preparation.py)

$ python gene_transfer.py wdir strain FastaPath

### TE transfer 

$ python TE_transfer.py wdir strain RefrencePath TEAnnotation
