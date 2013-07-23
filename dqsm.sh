#!/bin/sh


executable=./src/$1/dqsm
datadir=./data/legacy
data="-DeltaB $datadir/deltab.bin"
mask="-mask $datadir/mask.bin"
modelmap="-modelmap $datadir/models.bin"
maxiters="-maxiters $3"
model="-model m"
alpha="-alpha 0.5"
np="-np $2"
out="-out $outdir"

#previter="-x out41736/x_iter000925.bin"

today=$(date "+%d%m%Y")
outdir=out/$1/$today/
if [ -d $outdir ]; then
    rm -rf $outdir
else
    mkdir -p $outdir
fi
mkdir $outdir


#echo mpirun $np ./chimap $model $data $mask $modelmap $out $maxiters $alpha
mpirun $np $executable $model $data $mask $modelmap $out $maxiters $alpha
