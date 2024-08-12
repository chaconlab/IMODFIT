#!/bin/bash
echo "COMPILING & LINKING iMODFIT (Release v1.06)"
export IMOD_FLAGS=""
echo
echo "Entering: src"
cd src

# COMPILING LIBRARIES AND EXECUTABLES
for i in libtools libpdb libvolume  libnma libnmafit libgausscorr nmafit pdb2vol rmsd voltool 
 do
  echo "  Entering: $i"
  cd $i
   for j in *
   do
    if test -d $j
    then
     echo "    Entering: $j"
     cd $j
      echo Cleaning $i $j
      make clean
      echo Making $i $j
      make all
     cd ..
     echo "    Exiting: $j"
    fi
   done
  cd ..
  echo "  Exiting: $i"
 done

cd ..
echo "Exiting: src"
