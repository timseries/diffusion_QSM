dqsm (previously chimap) usage information.
===========================================

usage: ./dqsm <parameters>

Mandatory parameters:
---------------------

-DeltaB <deltab.bin>   : deltab data file
-mask <mask.bin>       : mask data file
-modelmap <models.bin> : model map data file
-out                   : output directory (must already exist)

Optional parameters:
--------------------

-model <s|m>                   : use (s)pherical model only, or (m)ixed models (ie 
                                 spherical and cylindrical). Default = m.
-chi <chi.bin>                 : initial chi array estimate. Default is all zeros.
-threshold <kernel threshold>  : minimum threshold for kernels. Any value below this 
                                 will be treated as zero. Default = 1e-6.
-mt <FA threshold>             : model FA threshold. A voxel with a value above this 
                                 will be modelled as a cylinder. Default = 0.2.
-rth <relative threshold>      : stopping criteria. Iterations will cease when the ratio 
                                 of the RMS of the change in solution to the RMS of the
                                 solution changes by less than this threshold. Default =
                                 1e-6.
-ath <absolute threshold>      : stopping criteria. Iterations will cease when the RMS of
                                 the change in the solution is less than this threshold.
                                 Default = 1e-16.
-maxiters <maximum iterations> : maximum number of iterations. Default = 1000.
-tau <tau>                     : The tau parameter in the dQSM cost function. Default =      
                                 0.15
-alpha <alpha>                 : The kappa parameter in the dQSM cost function. Default =
                                 0.75.
-x <x_iter*.bin>               : intermediate save state from previous runs (eg 
                                 x_iter000001.bin). Use this to continue calculations when 
                                 previous runs were prematurely stopped due to time outs.

Examples (on Avoca):
--------------------

Example 1: minimum parameters required

srun --ntasks-per-node=16 ./dqsm -DeltaB deltab.bin -mask mask.bin -modelmap models.bin -out output

Example 2: specifiying optional parameters

srun --ntasks-per-node=16 ./dqsm -DeltaB deltab.bin -mask mask.bin -modelmap models.bin -out output -model m -chi chi.bin -threshold 1e-6 -mt 0.2 -rth 1e-6 -ath 1e-16 -maxiters 1000 -tau 0.15 -alpha 0.75

Example 3: continuing a previous run that stopped prematurely

srun --ntasks-per-node=16 ./dqsm -DeltaB deltab.bin -mask mask.bin -modelmap models.bin -out output -x x_iter000900.bin

Outputs:
--------

All output is saved to the designated output directory (parameter -out). The outputs are:

chi.bin          : final susceptibility map
out.bin, out.m   : files for loading output into matlab
x_iter000000.bin : intermediate saved states during computation
x_iter000001.bin
x_iter...