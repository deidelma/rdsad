# rdsad

## Introduction

`rdsad` is a project aimed at simplifying the task of converting .rds files from [R](https://www.r-project.org) to [AnnData](https://github.com/scverse/anndata) file format suitable for use with [scanpy](https://github.com/scverse/scanpy).

## What is the problem this project aims to solve

This project arose from issues encountered in the [Ding Lab](https://junding.lab.mcgill.ca), which primarily uses Python to develop tooling for the analysis of single-cell RNA seq data.  However, much of the publicly available single-cell data is in rds format, produced using the excellent [Seurat](https://satijalab.org/seurat/) library. Understandably, when loaded onto public websites, this data is often encoded in `rds` files. Unfortunately, our experience has been that installing Seurat and the required libaries to convert the RDS files to `h5ad` format is fraught with problems associated with incompatabilities in library versions.

It turns out that this problem can be circumvented by carefully building a [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) environment with compatible libraries.  To ensure that all components be derived from a congruent set of libraries, the best approach is to create a clean Conda environment using the conda command line tool or equivalent (e.g. [mamba](https://github.com/mamba-org/mamba) or [pixi](https://github.com/prefix-dev/pixi) and then using a package source that is properly curated.

To accomplish this, it is essential to establish `conda-forge` as the channel to be used for downloading packages.  Assuming that "singlecell" is the environment name (use any name you choose), use the following steps.

1) create the environment:
```bash
conda create -n singlecell -c conda-forge
```
Once this is created, it is necessary to add a series of packages to the singlecell environment.

2) activate the environment:
```bash
conda activate singlecell
```
3) install mamba:

This is optional, as it is possible to use the `conda` command to install packages, but `mamba` runs much more quickly.

```bash
conda install mamba
```
This has the side effect of installing the current Python version, 3.12.

4) install scanpy:

```bash
mamba install scanpy
```
5) install R and needed R packages:

```bash
mamba install R, r-seurat, r-devtools, r-hdf5r
```
6) use R to load needed packages from github:

These replace the versions on conda-forge, which currently dod not work.

Open the R application by typing R at the command line.
```bash
R
```

Within R, enter the following commands:

```R
Library(Seurat)
Library(devtools)
devtools::install_github("satijalab/seurat-data")
devtools::install_github("mojaveazure/seurat-disk")
Library(SeuratData)
Library(SeuratDisk)
```
For each of the devtools commands, when prompted to update packages, press Return to accept the default.

7) if everything went well, you should now be able to use R to load and convert an `.rds` file into an `.h5ad` file. 

The chosen approach is based on the method outlined by Z.Q. Fang on this [website](https://zqfang.github.io/2020-04-28-seurat2scanpy/).  Although the author recommends MuDataSeurat to do the conversion, this does not seem to work with the latest libraries found on [conda-forge](https://conda-forge.org)).  Instead, it is recommended to use the approach proposed by the Satija Lab, which has authored Seurat.

First convert factor to character:
```R
# load the data from disk
data <- readRDS("myfile.rds")

# if desired, you can use DietSeurat to slim down the final result.  
# See https://zqfang.github.io/2020-04-28-seurat2scanpy/ for more information

# convert factor to character to avoid it being numbers in the h5ad file
i <- sapply(data@meta.data, is.factor)
data@meta.data[i] <- lapply(data@meta.data[i], as.chacter)

# convert the file
SaveH5Seurat(data, filename="myfile.h5seurat", overwrite=TRUE)
Convert("myfile.h5seurat", dest='h5ad', assay='RNA', overwrite=TRUE)
```
