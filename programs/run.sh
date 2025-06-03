#!/bin/bash -l

# Print current date
date

# Print name of node

ulimit -S -n 1000000

npart=( 2 4 8 16 24 32 ); 
output=./outputs/;
partsD=parts_${npart};
testd=~/bigdata/gkaur007/datasets/inputs/${partsD}/;  ##/scratch/gkaur007/inputs/${partsD}/;     ###~/bigdata/gkaur007/datasets/input
numLines=3072626;
input=~/bigdata/gkaur007/datasets/inputs/ok;   ##/scratch/gkaur007/inputs/adj_orkut.graph;
ftype=adjlist;
outpre=testFO;
outf=intGO32;
bsize=( 3072626 ); ## 3072626 2304470 1536313 768158 ); 
hideg=1000;
nthreads=32;
k=20;

for batchsize in ${bsize[@]};
do
   for nparts in ${npart[@]} ;
   do
        echo -e "1/4                               " >> ${output}${outf} ; 
        echo -e "\n*************************\n\n" >> ${output}${outf} ; 
        echo " /usr/bin/time -f "%P %M" ./kparts.bin ${input} ${ftype} ${numLines} ${hideg} ${nthreads} ${nparts} ${batchsize} ${k} ${testd}${outpre}${batchsize}${nparts}" >> ${output}${outf}; 
        ulimit -n 32768;
        srun ./kparts.bin ${input} ${ftype} ${numLines} ${hideg} ${nthreads} ${nparts} ${batchsize} ${k} ${testd}${outpre}${batchsize}${nparts} >> ${output}${outf} 2>&1; 
   done
done
