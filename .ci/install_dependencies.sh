#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  build-essential \
  cmake \
  git \
  libbz2-dev \
  libcurl4-gnutls-dev \
  liblzma-dev \
  libssl-dev \
  python3-pip \
  python3-setuptools \
  tabix \
  zlib1g-dev



if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root


#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
cp -s bcftools-1.10.2/bcftools .

#________________________ vt __________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .


#______________________vcflib _________________________________#
cd $install_root
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
make


#______________________ python ________________________________#
pip3 install tox "six>=1.14.0"
