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


### Do Init stage ###
cp TEMP_init CONTROL
echo "Running the simulation, Stage Init"
DLPOLY.Z
until [ $(grep -s "terminating" OUTPUT | wc | awk '{print $1}') -eq 1 ]
	do
		sleep 60
	done
cp REVCON CONFIG
mv OUTPUT OUTPUT_init
mv RDFDAT RDFDAT_init
mv REVCON REVCON_init
mv HISTORY HISTORY_init
#rm output_res HISTORY
echo '========================= DONE STAGE Init ========================='

### Do V1 stage ###
cp TEMP_V1 CONTROL
echo "Running the simulation, Stage V1"
DLPOLY.Z
until [ $(grep -s "terminating" OUTPUT | wc | awk '{print $1}') -eq 1 ]
	do
		sleep 60
	done
cp REVCON CONFIG
mv OUTPUT OUTPUT_V1
mv RDFDAT RDFDAT_V1
mv REVCON REVCON_V1
mv HISTORY HISTORY_V1
#rm output_res HISTORY
echo '========================= DONE STAGE V1 ========================='


### Do V4 stage ###
cp TEMP_V4 CONTROL
echo "Running the simulation, Stage V4"
DLPOLY.Z
until [ $(grep -s "terminating" OUTPUT | wc | awk '{print $1}') -eq 1 ]
        do
                sleep 60
        done
cp REVCON CONFIG
mv OUTPUT OUTPUT_V4
mv RDFDAT RDFDAT_V4
mv REVCON REVCON_V4
mv HISTORY HISTORY_V4
#rm output_res HISTORY
echo '========================= DONE STAGE V4 ========================='

exit 0

