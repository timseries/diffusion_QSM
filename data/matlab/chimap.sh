#!/bin/sh

mkdir output
mpirun -n 2 ./chimap -model m -DeltaB data.bin -mask mask.bin -modelmap models.bin -out output -maxiters 150