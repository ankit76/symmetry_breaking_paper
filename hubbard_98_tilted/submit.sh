#!/bin/bash
#SBATCH -A ucb-summit-sha
##SBATCH --qos ucb-summit-sha
##SBATCH --qos blanca-sha
#SBATCH --job-name hubbard98
#SBATCH --nodes 2
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --export=NONE


python sweepU.py > sweepU.out
