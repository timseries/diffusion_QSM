#!/bin/sh

#usage ./dqsm build_type number_of_nodes number_of_processes number_of_processes_per_node number_of_threads number_of_iterations
#build_type can be either "hybrid" or "legacy"
#example ./dqsm.sh "hybrid" 2 2 1 2 4 runs the executable in the "hybrid directory"
# with 2 nodes, 2 processes, 1 process per node, 2 threads, and 4 iterations

#environment variables
export HPM_OUTPUT_PROCESS=ROOT
export HPM_SCOPE=node
export HPM_ASC_OUTPUT=yes
export HPM_VIZ_OUTPUT=no
export HPM_METRICS=yes
export OUTPUT_ALL_RANKS=yes
export OMP_NUM_THREADS=$5

#previter="-x out41736/x_iter000925.bin"
today=$(date "+%d%m%Y")
outdir=out/$1/$today/n$2p$3ppn$4t$5i$6_all_profile
if [ -d $outdir ]; then
    echo "composite profile data directory exists, removing..."
    rm -rf $outdir
else
    mkdir -p $outdir
fi
for profile in mpi pomp hpm gmon
do
if [ $7 -eq 1 ]; then
make clean hybrid=1
make hybrid=1 omp=1 bluegene=1 debug=1 "$profile"_profile=1
#make hybrid=1 omp=1 debug=1
./dqsm.sh $1 $2 $3 $4 $5 $6 "$profile"_profile
else
#combine mode
cp out/$1/$today/n$2p$3ppn$4t$5i$6_"$profile"_profile/"$profile"* $outdir 
fi
done
make clean hybrid=1