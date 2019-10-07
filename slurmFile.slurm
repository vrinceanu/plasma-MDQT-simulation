#!/bin/bash
#SBATCH --job-name=testMD
#SBATCH --array=1-15
#SBATCH --partition=commons
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00

echo "My job ran on:"
echo $SLURM_NODELIST
if [[ -d $SHARED_SCRATCH/$USER && -w $SHARED_SCRATCH/$USER ]]
then
    srun runFile $SLURM_ARRAY_TASK_ID
fi
