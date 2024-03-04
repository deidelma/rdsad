FROM ubuntu:jammy  


RUN apt update 
RUN apt upgrade -y
RUN apt install curl -y
RUN apt install vim -y

RUN useradd -u 1234 meakins
USER meakins

WORKDIR /home/meakins

COPY ./rdsad-install.sh /home/meakins
COPY ./rdsad-install.R /home/meakins 
COPY ./install-conda.sh /home/meakins
COPY ./.gitignore /home/meakins
COPY ./LICENSE.txt /home/meakins
COPY ./README.md /home/meakins

USER root
RUN chmod +x rdsad-install.sh 
RUN chmod +x rdsad-install.R
RUN chmod +x install-conda.sh

RUN bash ./install-conda.sh 
ENV PATH "$PATH:/home/meakins/miniconda3/bin:/home/meakins/Miniconda3-latest-Linux-aarch64"

USER meakins 
RUN bash ./Miniconda3*.sh -b
COPY ./scripts /home/meakins/scripts
RUN conda init
RUN ./rdsad-install.sh 