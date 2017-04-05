#! /bin/sh

rm -rf espresso-5.1_modules_for_BandUP
tar -xvzf espresso-5.1_modules_for_BandUP.tgz
cd espresso-5.1_modules_for_BandUP

./configure FC=${1:-ifort} CC=${2:-icc} ${3:---enable-openmp} ${4:---disable-parallel}

make -j --max-load 2.5 bandup
