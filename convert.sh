#!/usr/bin/bash

#
# convert.sh
#
# script to coordinate the conversion of .rds files to .h5ad
#
# Copyright (c) 2024 David Eidelman. MIT License.
#

OVERWRITE=false

Help(){
    # echo "A script to convert Seurat .rds files into AnnData .h5ad files."
    echo "Usage: convert [-o|-h] filename"
    echo    filename needs to have the suffix .rds
    echo
    echo "options:"
    echo "o     Overwrite an existing .h5ad file."
    echo "h     Print this help."
    # echo "V     Print software version and license information."
    echo
}


#
# ensure filename is in the form *.rds
#
valid_filename(){
    [[ $1 =~ ^.*\.rds$ ]] 
}

# echo "Welcome to the rdsad file conversion program."
# echo

while getopts ":ho:" option; do
    case $option in
    h) # display help
        Help
        exit;;
    o) # allow overwrite
        echo "Overwriting output"
        OVERWRITE=true
        ;;
    \?) # invalid option
        echo "Invalid option"
        echo
        Help
        exit;;
    esac
done

FILENAME=${@: -1}
if [ -z $FILENAME ]; then
    # echo -n "Enter a filename:"
    # read -r FILENAME
    Help;
    exit 1
fi
if ! valid_filename $FILENAME; then
    echo Filename $FILENAME  is not valid.
    echo The Filename must end in .rds
    exit 1
fi
if ! [[ -e $FILENAME ]]; then
    echo Unable to find $FILENAME.  Did you enter the filename correctly?
    exit 1
fi
RSCRIPT="rdsad-convert.R"
if ! [[ -e $RSCRIPT ]]; then
    echo "Unable to find R script file $RSCRIPT. Did you install the rdsad system correctly?"
    exit 1
fi
if ! [[ -e ${FILENAME/.rds/.h5ad} ]] || $OVERWRITE; then
    echo "Converting the rds file."
    Rscript rdsad-convert.R "$FILENAME"
else 
    echo "Conversion already completed."
fi
INTERFILE=${FILENAME/".rds"/".h5seurat"}
if [[ -e $INTERFILE ]]; then
    echo "Removing intermediate file: $INTERFILE"
    rm $INTERFILE
fi

echo
echo rdsad has completed its tasks.