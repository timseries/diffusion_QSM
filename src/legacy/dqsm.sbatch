#!/bin/sh
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --job-name="dQSM"
#SBATCH --output="out%j.txt"

outdir=out$SLURM_JOB_ID
mkdir $outdir

datadir=../../data/legacy
data="-DeltaB $datadir/deltab.bin"
mask="-mask $datadir/mask.bin"
modelmap="-modelmap $datadir/models.bin"
out="-out $outdir"
maxiters="-maxiters 2"
model="-model m"
alpha="-alpha 0.5"
#previter="-x out41736/x_iter000925.bin"

time srun --ntasks-per-node=16 ./chimap $model $data $mask $modelmap $out $maxiters $previter $alpha
