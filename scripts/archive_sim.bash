#!/bin/bash

set -u
set -e
#set -v

simname=$1
archname=$2


echo "Making directory /nfs/data/ccsi/b0u/Kougarok/$archname"

mkdir -p /nfs/data/ccsi/b0u/Kougarok/$archname
echo "ncrcat $MYSCRATCH/$simname/run/*.h0.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h0.nc"
ncrcat $MYSCRATCH/$simname/run/*.h0.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h0.nc
echo "ncrcat $MYSCRATCH/$simname/run/*.h1.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h1.nc"
ncrcat $MYSCRATCH/$simname/run/*.h1.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h1.nc
echo "ncrcat $MYSCRATCH/$simname/run/*.h2.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h2.nc"
ncrcat $MYSCRATCH/$simname/run/*.h2.*.nc /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h2.nc

$HOME/plotting_py3_env/bin/python plot_kougarok.py /nfs/data/ccsi/b0u/Kougarok/$archname/${archname}_h2.nc --biomass --save_pftfile --save_figs /nfs/data/ccsi/b0u/Kougarok/$archname


