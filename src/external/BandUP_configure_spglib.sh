#! /bin/sh

rm -rf spglib-1.5.2
tar -xvzf spglib-1.5.2.tar.gz
cd spglib-1.5.2

export FC=ifort
export CC=icc
export CFLAGS='-openmp'
./configure FC=$FC CC=$CC

make
