#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

R CMD BATCH --no-save 20190730_time-memory.R 20190730_time-memory.Rout
