#!/bin/bash


CMD="./kparts.bin ${inputPath} ${fileFormattype} ${numVertices} ${hiDeg} ${nthreads} ${numPartitions} ${batchsize} ${klAlgoWindow} ${outputfilePath} >> ${output}${outf} 2>&1"; 
#-- CMD="./kpart --prodsm-compthreads=$1 --prodsm-vertexlist=$2 --prodsm-edgelist=$3 --prodsm-vertexmap=$4 --prodsm-iters=$5";

echo "Running: ${CMD}";
${CMD};
