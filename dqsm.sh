#!/bin/sh

#usage ./dqsm build_type number_of_nodes number_of_processes number_of_processes_per_node number_of_threads number_of_iterations
#build_type can be either "hybrid" or "legacy"
#example ./dqsm.sh "hybrid" 2 2 1 2 4 runs the executable in the "hybrid directory"
# with 2 nodes, 2 processes, 1 process per node, 2 threads, and 4 iterations
executable=./src/$1/dqsm
datadir=./data/legacy
data="-DeltaB $datadir/deltab.bin"
mask="-mask $datadir/mask.bin"
modelmap="-modelmap $datadir/models.bin"
maxiters="-maxiters $6"
model="-model m"
alpha="-alpha 0.5"
np="-np $3"
bluegenedir=/bgsys

#environment variables

export OMP_NUM_THREADS=$5

#previter="-x out41736/x_iter000925.bin"

today=$(date "+%d%m%Y")
outdir=out/$1/$today/n$2p$3ppn$4t$5i$6_$7
if [ -d $outdir ]; then
    echo "outdir exists, removing old data"
    rm -rf $outdir
else
    mkdir -p $outdir
fi
mkdir $outdir
out="-out $outdir"

#cross platform support
echo $bluegenedir
if [ -d $bluegenedir ]; then
    echo "we're making/running on bluegene"
#copy the sbatch template to the output directory
    cp ./dqsm.sbatch.template $outdir/dqsm.sbatch
#modify the sbatch file
    echo "#SBATCH --nodes="$2 >> $outdir/dqsm.sbatch
    echo "" >> $outdir/dqsm.sbatch
    echo "mkdir output" >> $outdir/dqsm.sbatch
    echo "executable=\"./dqsm\"" >> $outdir/dqsm.sbatch
    echo "data=\"-DeltaB ../../../.$datadir/deltab.bin\"" >> $outdir/dqsm.sbatch
    echo "mask=\"-mask ../../../.$datadir/mask.bin\"" >> $outdir/dqsm.sbatch
    echo "modelmap=\"-modelmap ../../../.$datadir/models.bin\"" >> $outdir/dqsm.sbatch
    echo "out=\"-out ./output\"" >> $outdir/dqsm.sbatch
    echo "maxiters=\"$maxiters\"" >> $outdir/dqsm.sbatch
    echo "model=\"$model\"" >> $outdir/dqsm.sbatch
    echo "alpha=\"$alpha\"" >> $outdir/dqsm.sbatch
    echo "" >> $outdir/dqsm.sbatch
    echo "time srun --ntasks-per-node="$4" \$executable \$model \$data \$mask \$modelmap \$out \$maxiters \$previter \$alpha" >> $outdir/dqsm.sbatch
#copy the executable to the output directory
    cp $executable $outdir
#change to that directory
    cd $outdir
#sbatch call
    sbatch dqsm.sbatch
else
#run locally
    echo "mpirun $np $executable $model $data $mask $modelmap $out $maxiters $alpha"
    mpirun $np $executable $model $data $mask $modelmap $out $maxiters $alpha
fi
