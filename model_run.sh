#!/bin/bash

# This is the script which controls the simulated annealing routine.
# It takes one parameter, -o (either true or false) indicating if
# you want to run the optimization or just the model using existing paramters.

# model setup directory
script_dir=`pwd`

usage() {                                      
  echo "Usage: $0 [ -o optimize: true or false ]" 1>&2
  exit 1 
}

while getopts ":o:" options; do
  case "${options}" in
    o) optimize=${OPTARG};;
    *) usage;;
  esac
done

# remove old files if they exist
files="phenograss *.mod *.txt "
for i in $files;
do
if [ -a $i ]; then
	rm $i
fi
done

# fresh compile code
gfortran -ffree-form \
 -ffree-line-length-200 \
 -g func.f90 inc.f90 phenograss.f90 sann.f90 main.f90 \
 -o phenograss

echo $optimize

# select routine
if [ "$optimize" == "true" ];then
  # run the optimization routine
  # full parameter output is stored in
  # summarizing_parameters.txt in the main
  # directory, either pick the last ones or
  # do weighted scheme 
  ./phenograss ./parameters/sites.txt o
else
  # run the model
  ./phenograss ./parameters/sites.txt r
fi

