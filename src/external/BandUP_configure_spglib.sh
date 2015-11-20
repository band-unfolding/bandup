#! /bin/sh

rm -rf spglib-1.5.2
tar -xvzf spglib-1.5.2.tar.gz
cd spglib-1.5.2

./configure FC=${1:-ifort} CC=${2:-icc} CFLAGS=${3:-''}

make
