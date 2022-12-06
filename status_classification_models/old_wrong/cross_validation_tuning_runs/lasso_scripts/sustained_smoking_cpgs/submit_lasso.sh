#!/bin/bash

#SBATCH -p exacloud                     # partition (queue)
#SBATCH -N 1                            # number of nodes
#SBATCH -n 16                           # number of cores
#SBATCH --mem 250000                    # memory pool for all cores
#SBATCH -t 0-24:00                      # time (D-HH:MM)
#SBATCH -o lasso_%A_%a.out              # Standard output
#SBATCH -e lasso_%A_%a.err              # Standard error
#SBATCH --array=1-200                   # sets number of jobs in array

# activate conda environment
source activate tidymodels

echo "SLURM_JOBID: " $SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID

# create hyperparameters args file
#Rscript make_hyperparams_lasso.R
# arguments file
dir=$(pwd)
args_file=$dir/lasso_hyperparam_args_file

# set variable "ARGS" to be the output of the correct line from args_file (matches SLURM_ARRAY_TASK_ID)
ARGS=`head -$SLURM_ARRAY_TASK_ID $args_file | tail -1`
echo $ARGS
# set variables based on the output of the line ARGS
IFS=: read H1 H2 SEED <<< $ARGS

echo "H1: " $H1
echo "H2: " $H2
echo "SEED: " $SEED

# execute the modeling script on each hyperparameter combo
Rscript run_lasso.R -a $H1 -b $H2 -s $SEED

## manualyl check for failed jobs
#grep -i "^error" *.err | cut -d : -f1
## copy those lines from args file to a new args file
#sed -n '27p;69p;121p;159p' lasso_hyperparam_args_file > lasso_hyperparam_args_file_v2
## rerun the new jobs by changing the "--array" number and the args file name


#wait

# Check for any failed jobs by grepping for "error" in the .err files
#num_errs=$(grep -i "^error" *.err | wc -l)
#if num_errs > 0
#   #grep -i "^error" *.err | cut -d : -f1 | sed 's/.err//g' | sed  "s/$SLURM_ARRAY_JOB_ID//g" | sed 's/lasso__//g'
#   grep -i "^error" *.err | cut -d : -f1 > failed_jobs.txt
#fi
# summarize performance metrics from all CV runs
#Rscript summarize_performance_metrics_lasso.R

# create best performing hyperparams file and print best performing hyperparams
