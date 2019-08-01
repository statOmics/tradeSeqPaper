#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1

R CMD BATCH --no-save generate_dyntoy_dataset.R generate_dyntoy_dataset.Rout
