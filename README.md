# MODCAN: driver gene identification based on multi-omics features and differential co-association networks for tumor subtypes

### This is the original repository for the MODCAN paper. 

**Requirements**

Python 3.13

```
h5py>=3.14.0
numpy>=2.3.1
networkx>=3.5
pandas>=2.3.0
scikit-learn>=1.7.0
scipy>=1.16.0
torch>=2.7.1
torch-geometric>=2.6.1
```

# **Input**

## Protein Protein Interaction (PPI) Network

MODCAN utilizes the PPI network, STRINGv12, as the underlying network. You can download this network file from https://cn.string-db.org/cgi/download.

The file should be stored at ./data/network/string_full_v12.txt.

Example:
```
protein_1	protein_2	score
ARF5	        CYTH2	        471
ARF5	        GGA1	        594
ARF5	        GOSR2	        303
ARF5	        MET	        347
...
```

## Genomic Data
Gene expression, promoter methylation, and copy number variant (CNV) data for cancer patients, sourced from TCGA, can be retrieved through the UCSC Xena Browser: https://xenabrowser.net/datapages/.

Single nucleotide variation (SNV) data can be downloaded via the R package "TCGAbiolinks".

Example:
```
library('TCGAbiolinks')

projects <- sort(getGDCprojects()$project_id)
cancer_type <- "TCGA-****"    # "TCGA-BLCA" "TCGA-BRCA" "TCGA-COAD" "TCGA-HNSC" "TCGA-KIRC" "TCGA-KIRP" "TCGA-LUAD" "TCGA-LUSC" "TCGA-PRAD" "TCGA-UCEC"  
Project_Summary <- getProjectSummary(cancer_type)

# download SNV data
query <- GDCquery(
  project = cancer_type,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  access = "open",
  legacy = FALSE
)
GDCdownload(query)

save_path <- "D:/MODCAN/TCGA_data/IHGC/TCGA_hg38_SNV/"
# dir.create(save_path)
GDCprepare(query, save = TRUE, save.filename = paste0(save_path, cancer_type, "_hg38_SNV.Rdata"))
```

### All tab-delimited genomic data should be pre-processed and stored in the following format. For a given cancer type, all genomic data for tumor samples should share the same row and column names.

## SNV Data
Store the SNV data in the folder:  ./data/TCGA_hg38_SNV/snv_matrix/, including a gene ID file, a sample ID file and an SNV matrix file.

Example:
```
1. Gene ID File (e.g., TCGA-BLCA_snv_gene.txt): Contains gene IDs that correspond to the rows of the SNV matrix.
...
EPHA2
USP24
SERBP1
GJA5
...

2. Sample ID File (e.g., TCGA-BLCA_snv_sample.txt): Contains sample IDs that correspond to the columns of the SNV matrix.
...
TCGA-XF-AAMH-01A
TCGA-GD-A2C5-01A
TCGA-DK-AA6W-01A
TCGA-FD-A6TF-01A
...

3. SNV Matrix File (e.g., TCGA-BLCA_snv_matrix.txt): Contains the SNV data matrix.
...
1	0	0	1	1	...
0	0	1	0	1	...
1	1	0	0	0	...
0	1	0	1	0	...
...
```

### Store the gene expression, promoter methylation, and CNV data in separate folders. These should be formatted similarly to the SNV data.

Gene expression:	./data/TCGA_UCSC_EXP/TCGA-BLCA/

Promoter methylation: ./data/TCGA_UCSC_MET/TCGA-BLCA/

CNV: ./data/TCGA_UCSC_CNV/cnv_matrix/

## Gene Expression and Promoter Methylation Data for Normal Samples

For TCGA-BLCA, the data should be stored in the folder: ./data/TCGA_UCSC_normal/TCGA-BLCA/, including  a gene ID file, a sample ID file and two data matrices.

These files should be in the same format as the SNV data.

## Clinical Data
Clinical data can also be downloaded from the UCSC Xena Browser: https://xenabrowser.net/datapages/.

For TCGA-BLCA, the clinical data should be stored in ./data/survival/TCGA-BLCA.survival.tsv.

Example:
```
sample	                OS	_PATIENT        OS.time
TCGA-E7-A8O8-01A	0	TCGA-E7-A8O8	13
TCGA-GC-A4ZW-01A	0	TCGA-GC-A4ZW	15
TCGA-E7-A5KE-01A	0	TCGA-E7-A5KE	17
TCGA-4Z-AA80-01A	1	TCGA-4Z-AA80	19
...
```

## Reference Driver Genes

The reference driver genes file is located at: ./data/reference/CGC.txt.

# **Run**

Example:

```
1. cd ./R_code/
	Run ./Cluster.R

2. cd ./src/
	python train_Mgcn.py -e 1000 -lr 0.001 -hd [64, 128] -lm 2 -wd 5e-4 -do 0.5 -wt 2 -ns 0.2 -d '../h5_file/TCGA-BLCA_string_full_v12_0.7.h5' -cv 10 -dr 'cgc'
```
