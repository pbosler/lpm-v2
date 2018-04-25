#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --account=fy150039
#SBATCH --job-name=gHills-lpmI
#SBATCH --partition=ec
#SBATCH --mail-type=ALL

module purge
module load intel/15.0
module load gnu/4.9.2
module load openmpi-intel/1.8

LPM_ROOT=$HOME/lpm-v2
#OUTPUT_LOC=$HOME/modelData/transport
OUTPUT_LOC=/fscratch/pabosle/modelData/transport/timingIndirect

nnodes=$SLURM_JOB_NUM_NODES
# nnodes=1
cores=8
nnp=$(($cores*$nnodes))

JOB_FILENAME=gHillsIndOutput_np${nnp}_parRM.txt

cp $LPM_ROOT/build/advectGaussHillsSphere.exe $OUTPUT_LOC/.
cd $OUTPUT_LOC

rm -f $JOB_FILENAME

for i in `seq 1 7`
do
cat <<EOF > gHillsInd.namelist
&meshDefine
	faceKind = 3 
	initNest = ${i}
	amrLimit = 0
/case 

&timestepping
	dt = 0.0125
	tfinal = 5.0
	remeshInterval = 20
	useDirectRemesh = .FALSE.
/

&fileIO
	outputDir = '${OUTPUT_LOC}'
	outputRoot = 'gHills_indirect_np${nnp}'
	frameOut = 400
/
EOF

mpirun --bind-to core --npernode $cores -np $nnp advectGaussHillsSphere.exe gHillsInd.namelist 2>&1 | tee -a $OUTPUT_LOC/$JOB_FILENAME

done

