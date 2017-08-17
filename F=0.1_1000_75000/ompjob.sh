#!/bin/csh

### Set the job name
#PBS -N F=0.1_1000_75000

### Request email when job begins and ends
#PBS -m bea

### Specify email address to use for notification.
#PBS -M huitang@email.arizona.edu

### Specify the PI group for this job
### List of PI groups available to each user can be found with "va" command
#PBS -W group_list=lmcguire

### Set the queue for this job as windfall or standard (adjust ### and #)
#PBS -q standard
###PBS -q windfall

### Set the number of cores (cpus) and memory that will be used for this job
### When specifying memory request slightly less than 6GB memory per ncpus for standard node
### Some memory needs to be reserved for the Linux system processes
#PBS -l select=1:ncpus=28:mem=10gb

### Specify "wallclock time" required for this job, hhh:mm:ss
#PBS -l walltime=48:00:00

### Specify total cpu time required for this job, hhh:mm:ss
### total cputime = walltime * ncpus
#PBS -l cput=1536:00:00

### Load required modules/libraries if needed (blas example)
### Use "module avail" command to list all available modules
### NOTE: /usr/share/Modules/init/csh -CAPITAL M in Modules
###source /usr/share/Modules/init/csh
###module load blas

### set directory for job execution, ~netid = home directory path
cd ~huitang/F=0.1_1000_75000

###
setenv OMP_NUM_THREADS 28
unlimit

### run your executable program with begin and end date and time output
date
./2dmuscl_hr_omp
date
