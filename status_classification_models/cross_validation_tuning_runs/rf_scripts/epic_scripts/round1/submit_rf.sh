#!/bin/bash

#SBATCH -p exacloud                     # partition (queue)
#SBATCH -N 1                            # number of nodes
#SBATCH -n 16                           # number of cores
#SBATCH --mem 250000                    # memory pool for all cores
#SBATCH -t 0-24:00                      # time (D-HH:MM)
#SBATCH -o rf_%A_%a.out                 # Standard output
#SBATCH -e rf_%A_%a.err                 # Standard error
#SBATCH --array=1-100                   # sets number of jobs in array

# activate conda environment
source activate tidymodels

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# create hyperparameters args file
Rscript make_hyperparams_rf.R
# arguments file
dir=$(pwd)
args_file=$dir/rf_hyperparam_args_file

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`
echo $ARGS
# set variables based on the output of the line ARGS
IFS=: read H1 H2 SEED <<< $ARGS

echo "H1: " $H1
echo "H2: " $H2
echo "SEED: " $SEED

# execute the modeling script
Rscript run_rf.R -a $H1 -b $H2 -s $SEED
