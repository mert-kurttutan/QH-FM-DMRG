#!/bin/bash
#SBATCH --job-name FHH_SU2
#SBATCH --comment "FHH Model - QH Ferromagnetism - Convergence"
#SBATCH --time=3-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Mert.Kurttutan@physik.uni-muenchen.de
#SBATCH --chdir=.
#SBATCH --output=/project/th-scratch/m/Mert.Kurttutan/slurm/DMRG.%j.%N.out
#SBATCH --constraint=avx
#SBATCH --ntasks=4
#SBATCH --mem=128GB
#SBATCH --partition=th-cl,cluster,th-ws


# input: Lx Ly alpha U N S(spin) PBC job

source /project/theorie/s/Sam.Mardazad/Group/syten_inc.sh

#python ./run_DMRG.py $1 $2 $3 $4 $5 $6 $7
#python ./convergence_n_arr.py $1 $2 $3 $4 $5 $6 $7
#python ./convergence_cur_arr.py $1 $2 $3 $4 $5 $6 $7
#python ./sCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7
python ./nCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7
