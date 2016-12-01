#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:18:00
#SBATCH --account=fy150039
#SBATCH --job-name=gHills-lpmD
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

module purge
module load intel/15.0
module load gnu/4.9.2
module load openmpi-intel/1.8

nnodes=$SLURM_JOB_NUM_NODES
cores=8
nnp=$(($cores*$nnodes))

LPM_ROOT=$HOME/lpm-v2
#OUTPUT_LOC=$HOME/modelData/transport
OUTPUT_LOC=/fscratch/pabosle/modelData/transport/timingDirect
JOB_FILENAME=gHillsDirectOutput_np${nnp}_parRM.txt

cp $LPM_ROOT/build/advectGaussHillsSphere.exe $OUTPUT_LOC/.
cd $OUTPUT_LOC

rm -f $JOB_FILENAME

for i in `seq 1 7`
do
cat <<EOF > gHillsDirect.namelist
&meshDefine
	faceKind = 3 
	initNest = ${i}
	amrLimit = 0
/case 

&timestepping
	dt = 0.0125
	tfinal = 5.0
	remeshInterval = 20
	useDirectRemesh = .TRUE.
/

&fileIO
	outputDir = '${OUTPUT_LOC}'
	outputRoot = 'gHills_direct_np${nnp}'
	frameOut = 400
/
EOF

mpirun --bind-to core --npernode $cores -np $nnp advectGaussHillsSphere.exe gHillsDirect.namelist 2>&1 | tee -a $OUTPUT_LOC/$JOB_FILENAME

done

