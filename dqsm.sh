#!/bin/sh

if [ $1 = "hybrid" ]; then
    executable="src/hybrid/dqsm"
fi
outdir=out/$1
if [ -d $outdir ]; then
    rm -rf $outdir
fi

if [ -d out ]; then
    echo "out exists"
else
    mkdir out
fi

mkdir $outdir

executable=./src/$1/dqsm
datadir=./data/legacy
data="-DeltaB $datadir/deltab.bin"
mask="-mask $datadir/mask.bin"
modelmap="-modelmap $datadir/models.bin"
out="-out $outdir"
maxiters="-maxiters 2"
model="-model m"
alpha="-alpha 0.5"
np="-np 2"
#previter="-x out41736/x_iter000925.bin"

#echo mpirun $np ./chimap $model $data $mask $modelmap $out $maxiters $alpha
mpirun $np $executable $model $data $mask $modelmap $out $maxiters $alpha
