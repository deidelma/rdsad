#!/usr/bin/env bash

echo "rdsad-install"
echo ""
echo "This script installs the necessary packages for converting betweeen rds and h5ad files."
echo ""
echo "It is assumed that the script is excuting in a new, empty Conda environment that has"
echo "been created with the command:"
echo ""
echo "conda create -n myenv -c conda-forge"
echo ""
echo "'myenv' is an environment name of your choosing and conda-forge is the channel"
echo "that conda uses to retrieve packages."
echo ""
echo "Before executing this script, make sure to activate the new environment with the command:"
echo ""
echo "conda activate myenv"
echo ""
conda install -y mamba -c conda-forge
echo "mamba installed"
mamba install -y scanpy -c conda-forge
echo "scanpy installed"
mamba install -y -c conda-forge 'jupyterlab>=4.1.0,<5.0.0a0' jupyterlab-lsp
echo "jupyter installed"
mamba install -y -c conda-forge python-lsp-server r-languageserver
echo "jupyter language servers installed"
mamba install -y -c conda-forge rpy2
echo "rpy2 installed"

# mamba install -y R r-seurat r-devtools r-hdf5r -c conda-forge
conda install -y R -c conda-forge 
conda install -y r-seurat -c conda-forge 
conda install -y r-devtools -c conda-forge
# conda install -y r-hdf5 -c conda-forge
echo "R and associated packages installed"
./rdsad-install.R
echo "R packages installed from github"
echo ""
echo "Installation complete."
echo ""


