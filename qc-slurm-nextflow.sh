#!/bin/bash
#SBATCH -p WORK # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH -c 2 # number of cpus per task
#SBATCH --mem 10 # memory pool for all cores in Gb
#SBATCH -t 10:00:30 # max time (HH:MM:SS)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

echo "### START - the job is starting at"
date
starttime=`date +"%s"`
echo
echo "the job is running on the node $SLURM_NODELIST"
echo "job number $SLURM_JOB_ID"
echo "STAT:jobName:$SLURM_JOB_ID"
echo "STAT:exechosts:$SLURM_NODELIST"
echo

cd $PWD

#source /apps/conda3.2.sh
module load NextFlow/latest
nextflow -Dcapsule.log=verbose info
PROFILE=conda nextflow run ./main.nf -entry test_QC -c ./nextflow.config -resume

echo "####END job finished"
endtime=`date +"%s"`
duration=$((endtime - starttime))
echo "STAT:startTime:$starttime"
echo "STAT:doneTime:$endtime"
echo "STAT:runtime:$duration"
