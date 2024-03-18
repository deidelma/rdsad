#!/usr/bin/bash
#
# install-conda-interactive.sh
#
# Shell script to install conda interactively.
# Checks for active environment and handles reinstall correctly.
# Assumes that the environment name is 'rdsad'
#
# copyright (c) 2024 David Eidelman, MIT License.
#

if [ -z ${CONDA_PREFIX} ]; then
    echo "rdsad conda environment not detected." 
    read -p "Attempt to install enviornment? (Y/N): " confirm 
    if [[ $confirm == [yY] ]]; then
        conda create -n rdsad -c conda-forge python=3.11 -y
    else
        exit 1
    fi
else
    echo "rdsad conda environment active.  If you wish to reinstall, it must be deleted first."
    echo "Use 'conda env remove rdsad --all' to delete the current version."
    exit 0
fi

