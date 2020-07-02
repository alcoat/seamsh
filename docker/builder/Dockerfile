FROM ubuntu:20.04

env DEBIAN_FRONTEND=noninteractive 
#linux
RUN apt update && apt install -y git python3 make cmake gcc
#windows
RUN apt update && apt install -y mingw-w64 mingw-w64-tools
run echo "\
SET(CMAKE_SYSTEM_NAME Windows)\n\
SET(CMAKE_C_COMPILER x86_64-w64-mingw32-gcc)\n\
SET(CMAKE_CXX_COMPILER x86_64-w64-mingw32-g++)\n\
SET(CMAKE_RC_COMPILER x86_64-w64-mingw32-windres)\n\
SET(CMAKE_Fortran_COMPILER x86_64-w64-mingw32-gfortran)\n\
SET(CMAKE_FIND_ROOT_PATH /mingw64)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)\n\
SET(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)\n"\
> /cmake-mingw

RUN apt update && apt install -y git python3-setuptools python3-wheel twine


### osxcross

RUN apt-get update && apt-get install -y clang patch libssl-dev liblzma-dev libxml2-dev llvm-dev uuid-dev zlib1g-dev


RUN useradd -ms /bin/bash validator
RUN mkdir -p /osxcross && chown validator /osxcross 
USER validator
RUN git clone https://github.com/tpoechtrager/osxcross.git /osxcross
COPY MacOSX10.15.sdk.tar.xz /osxcross/tarballs/MacOSX10.15.sdk.tar.xz
RUN cd /osxcross && SDK_VERSION=10.15 UNATTENDED=1 ./build.sh
ENV PATH=/osxcross/target/bin/:$PATH
ENV OSXCROSS_HOST=x86_64-apple-darwin19

VOLUME ["/etc/gitlab-runner"]

WORKDIR /home/validator
