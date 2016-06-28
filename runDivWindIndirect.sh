#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:60:00
#SBATCH --account=fy150039
#SBATCH --job-name=adv-divWindGHills-lpmI
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

LPM_ROOT=$HOME/lpm-v2
#OUTPUT_LOC=$HOME/modelData/transport
OUTPUT_LOC=/fscratch/pabosle/modelData/transport
JOB_FILENAME=divWindGHillsOutputInd.txt

cd $LPM_ROOT/build

rm -f $JOB_FILENAME

for i in `seq 1 7`
do
cat <<EOF > divWindGHillsInd.namelist
&meshDefine
	faceKind = 3 
	initNest = ${i}
	amrLimit = 0
/

&timestepping
	dt = 0.0125
	tfinal = 5.0
	remeshInterval = 20
	useDirectRemesh = .FALSE.
/

&fileIO
	outputDir = '${OUTPUT_LOC}'
	outputRoot = 'divWind_gHills_indRemesh'
	frameOut = 100
/
EOF

mpirun -np 8 advectGHillsDivWind.exe divWindGHillsInd.namelist 2>&1 | tee -a $JOB_FILENAME

done

