#!/bin/sh
#SBATCH -n 1
#SBATCH -p cpu-6230R
#SBATCH --mail-type=BEGIN,END
#SBATCH -J PMF

#module load intel/19.1.3.304
#module load impi/2019.9.304
module load dl-poly/4.09_MahaGaro-intel



### All calculation parameters go here ###
CUWOD_SS4=$(pwd)
STEMP_SS4='4000'
NSTEPS_SS4='38'


### Calculation initialization ###
cd $CUWOD_SS4


DLPOLY.Z

cp REVCON CONFIG
mv OUTPUT OUTPUT_init
mv RDFDAT RDFDAT_init
mv REVCON REVCON_init
mv HISTORY HISTORY_init
#rm output_res HISTORY
echo '========================= DONE STAGE Init ========================='


DLPOLY.Z

mv OUTPUT OUTPUT_V4
mv RDFDAT RDFDAT_V4
mv REVCON REVCON_V4
mv HISTORY HISTORY_V4
#rm output_res HISTORY
echo '========================= DONE STAGE Init ========================='

exit 0


