FROM ubuntu:jammy  


RUN apt update 
RUN apt upgrade -y
RUN apt install curl -y
RUN apt install vim -y

RUN useradd -u 1234 meakins
USER meakins

WORKDIR /home/meakins

COPY . /home/meakins

RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod -v +x Miniconda3*.sh 
RUN bash ./Miniconda3*.sh -b

ENV PATH "$PATH:/home/meakins/miniconda3/bin"

USER meakins 

RUN conda init
RUN ./scripts/rdsad-install.sh 