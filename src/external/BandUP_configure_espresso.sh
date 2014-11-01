#! /bin/sh

rm -rf espresso-5.1
tar -xvzf espresso-5.1.tgz
cd espresso-5.1

export FC=ifort
export CC=icc
./configure FC=$FC CC=$CC

make bandup
