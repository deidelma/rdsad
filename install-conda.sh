#
# install_conda.sh
# 
# shell script to install conda on either Intel or Apple hardware
# designed to be called from a Docker container build or equivalent
# assumes that curl is already installed
#
# copyright (c) 2024 David Eidelman, MIT License
#

OS_NAME=$(uname -p);
echo "OS_NAME is $OS_NAME"

if [[ $OS_NAME == 'aarch64' ]]; then
    echo "Downloading Apple Silicon version of miniconda"
    curl -O "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-aarch64.sh"
else
    echo "Downloading x86_64 version of miniconda"
    curl -O "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
fi;
echo "Unpacking miniconda."
chmod -v +x Miniconda3*.sh 
# ./Miniconda3*.sh -b

