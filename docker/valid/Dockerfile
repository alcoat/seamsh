FROM ubuntu:20.04

RUN apt update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt update &&  apt install -y git python3 python3-numpy gmsh python3-scipy python3-petsc4py python3-pip gdal-bin

VOLUME ["/etc/gitlab-runner"]

RUN useradd -ms /bin/bash validator
USER validator
WORKDIR /home/validator
