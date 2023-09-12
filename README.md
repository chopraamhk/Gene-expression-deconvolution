# CIBERSORTX - Deconvolution Steps

STEP1:
Get the .hdf5 (anndata file - one of single cell file format) file from existing single cell study 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165824
  
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

```{r}
seurat_anndata@assays$RNA
```  

```{r}
here, the code to make the file right
```

How does the output_ref.txt looks like <tab delimited file>- 
```
Genes	VSMC_II	VSMC_I	VSMC_I	Endothelial_I
RP11-34P13.3	0	0	0	0
FAM138A	0	0	0	0
OR4F5	0	0	0	0
RP11-34P13.7	0	0	0	0
```
How does the .txt looks like <tab delimited file>- 
```
Genes	GTEX-111YS	GTEX-1122O	GTEX-1128S	GTEX-117XS
DDX11L1	0	0	0	0
WASH7P	3.652	3.402	7.733	3.057
MIR6859-1	0	0	0	0
MIR1302-2HG	0	0	0.0659	0
```
#Running cibersortx
```
sudo docker run -v /home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/results:/src/data -v /home/mchopra/Documents/PhD-Year1/deconvolution/Deconvolution_results/results:/src/outdir cibersortx/fractions --username m.chopra1@nuigalway.ie --token 32c3bd33e2bdbe49a89c001a948cf2e5 --single_cell TRUE --refsample output_ref.txt --mixture gtex_aorta.txt --perm 100
```
