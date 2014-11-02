#! /bin/sh

rm -rf espresso-5.1_modules_for_BandUP
tar -xvzf espresso-5.1_modules_for_BandUP.tgz
cd espresso-5.1_modules_for_BandUP

export FC=ifort
export CC=icc
./configure FC=$FC CC=$CC

make bandup
