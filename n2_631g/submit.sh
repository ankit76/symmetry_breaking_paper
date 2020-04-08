#!/bin/bash
#SBATCH -A ucb-summit-sha
##SBATCH --qos ucb-summit-sha
##SBATCH --qos blanca-sha
#SBATCH --job-name n2
#SBATCH --nodes 1
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --export=NONE


python pes.py > pes.out
