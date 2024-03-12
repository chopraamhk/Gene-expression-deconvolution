# CIBERSORTX - Deconvolution Steps
https://open.bioqueue.org/home/knowledge/showKnowledge/sig/cibersortx-fractions#--fraction

STEP1:
Get the .hdf5 (anndata file - one of single cell file format) file from existing single cell study 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165824

NOTE: MAKE SURE THAT YOU ARE CONSIDERING THE CORRECT SINGLE-CELL STUDY AS IT WILL BE THE REFERENCE
  
#load libraries
```{r}
library(biomaRt)
library(Seurat)
library(SeuratDisk)
library(tidyverse) #to manipulate the data
#library(HDF5Array)
#library(BiocManager)
#library(zellkonverter) #zellKonverter is a bioc package that can do the job as well
```

#setwd
```{r}
setwd("/home/mchopra/Documents/PhD-Year1/deconvolution/GSE165824_ascending_descending_human_aorta_v1")
```

#H5ad -> anndata object ##Data has been downloaded from NCBI Geo GSE165824.
#if want to load other types of data check -> https://www.youtube.com/watch?v=3xcTpqQzUwQ 
#Step 1: convert AnnData object to an h5Seurat file --> this will create a h5Seurat file

```{r}
#file <- readH5AD("GSE165824_ascending_descending_human_aorta_v1.h5ad", use_hdf5 = FALSE)  #if using zellKonverter
Convert("GSE165824_ascending_descending_human_aorta_v1.h5ad", dest = "h5seurat", overwrite = TRUE)
```

STEP2: Load h5Seurat file into a Seurat Object
```{r}
#file
seurat_anndata <- LoadH5Seurat("GSE165824_ascending_descending_human_aorta_v1.h5seurat")
#file@assays@data$X
```

```{r}
seurat_anndata
str(seurat_anndata)
#As it has two lists, we use RNA and saved it into variable called cts.
```

Step3: See the data
```{r}
seurat_anndata@assays$RNA
```  

```{r}
seurat_anndata@assays$RNA
table(seurat_anndata$donor_id)
table(seurat_anndata$biosample_id)
```

```{r}
aoA_samples <- subset(seurat_anndata, subset =  biosample_id %in% c("Ao4A", "Ao8A", "Ao12A"))
mtrx <- aoA_samples@assays$RNA@counts
```

```{r}
colnames(mtrx) <- aoA_samples$cell_type_leiden
```

```{r}
table(seurat_anndata$cell_type_leiden)
```

```{r}
colnames(mtrx)
```

```{r}
colnames(mtrx) <- gsub("00. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("03. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("04. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("05. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("01. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("02. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("09. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("08. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("06. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("10. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("13. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("07. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("11. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("12. ", "", colnames(mtrx))
colnames(mtrx) <- gsub("13. ?", "unknown", colnames(mtrx))
colnames(mtrx)
```

```{r}
write.csv(mtrx, "mtrx.csv", quote = FALSE, sep = "\t")
```

```{r}
output_ref <- read.csv("/home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/Sc_reference/mtrx.csv", check.names = FALSE)
```

```{r}
output_ref[1:5, 1:5]
```

```{r}
colnames(output_ref) <- sub("\\.", "_", colnames(output_ref))
colnames(output_ref) <- sub("\\..*", "", colnames(output_ref))
colnames(output_ref)[1] <- "Genes"
colnames(output_ref)
colnames(output_ref) <- sub("Pericyte_\\d+$", "Pericyte", colnames(output_ref))
colnames(output_ref) <- sub("Macrophage_\\d+$", "Macrophage", colnames(output_ref))
colnames(output_ref) <- sub("Lymphocyte_\\d+$", "Lymphocyte", colnames(output_ref))
colnames(output_ref) <- sub("Neuron_\\d+$", "Neuron", colnames(output_ref))
colnames(output_ref) <- sub("Mesothelial_\\d+$", "Mesothelial", colnames(output_ref))
column_names <- c(colnames(output_ref))

remove.trailing.numbers <- function(column_name) { gsub("\\.\\d+$", "", column_name)}

new_column_names <- sapply(column_names, remove.trailing.numbers)
print(new_column_names)
colnames(output_ref) <- new_column_names
output_ref[1:5, 1:5]
```

```{r}

colnames(output_ref) <- sub("\\..*", "", colnames(output_ref))
output_ref[1:5, 1:5]
```

```{r}
write.table(output_ref, "/home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/Sc_reference/output_ref.txt", sep = "\t", quote = FALSE, row.names = FALSE)
```

How does the output_ref.txt looks like - tab delimited file- 
```
Genes	VSMC_II	VSMC_I	VSMC_I	Endothelial_I
RP11-34P13.3	0	0	0	0
FAM138A	0	0	0	0
OR4F5	0	0	0	0
RP11-34P13.7	0	0	0	0
```
How does the .txt looks like tab delimited file- 
```
Genes	GTEX-111YS	GTEX-1122O	GTEX-1128S	GTEX-117XS
DDX11L1	0	0	0	0
WASH7P	3.652	3.402	7.733	3.057
MIR6859-1	0	0	0	0
MIR1302-2HG	0	0	0.0659	0
```
#Pull docker image 
```
docker pull cibersortx/fractions
```

#Running cibersortx
```
sudo docker run -v /home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/results:/src/data -v /home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/results:/src/outdir cibersortx/fractions --username <registered_email> --token <token from website> --single_cell TRUE --refsample output_ref.txt --mixture gtex_aorta.txt --perm 100
```
