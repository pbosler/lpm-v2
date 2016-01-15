#!/bin/bash

EXE_DIR=/Users/pabosle/lpm-v2/build
OUTPUT_DIR=/Users/pabosle/modelData/tests

pow=0.6
outFile=$OUTPUT_DIR/convTestResultsPow_$pow.txt

rm -f $outFile

cd $EXE_DIR

for (( i=3; i < 7; i++))
do
cat <<EOF > pseSphereConvergence.namelist
&meshDefine
	faceKind = 3
	initNest = $i
	amrLimit = 0
	radius = 1.0
	psePower = $pow
/
&fileIO
	outputDir='${OUTPUT_DIR}'
	outputRoot='spherePSEConvTest_pow${pow}'
/

EOF

mpirun -np 10 spherePSETest.exe pseSphereConvergence.namelist 2>&1 | tee -a $outFile

done

for (( i=3; i < 7; i++))
do
cat <<EOF > pseSphereConvergence.namelist
&meshDefine
	faceKind = 4
	initNest = $i
	amrLimit = 0
	radius = 1.0
	psePower = $pow
/
&fileIO
	outputDir='${OUTPUT_DIR}'
	outputRoot='spherePSEConvTest_pow${pow}'
/

EOF

mpirun -np 10 spherePSETest.exe pseSphereConvergence.namelist 2>&1 | tee -a $outFile

done

