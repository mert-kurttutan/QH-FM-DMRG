#!/bin/bash
#SBATCH --job-name FHH_SU2
#SBATCH --comment "FHH Model - QH Ferromagnetism"
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=Mert.Kurttutan@physik.uni-muenchen.de
#SBATCH --chdir=/home/m/Mert.Kurttutan/Academia/Codes/Physics/Projects/qh_fm_01/codes/excs
#SBATCH --output=/project/th-scratch/m/Mert.Kurttutan/slurm/DMRG.%j.%N.out
#SBATCH --constraint=avx
#SBATCH --ntasks=4
#SBATCH --mem=128GB
#SBATCH --partition=th-cl,cluster

__dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# input: Lx Ly alpha U N S(spin) PBC job

source /project/theorie/s/Sam.Mardazad/Group/syten_inc.sh

#python ${__dir}/run_DMRG.py $1 $2 $3 $4 $5 $6 $7
#python ${__dir}/convergence_n_arr.py $1 $2 $3 $4 $5 $6 $7
#python ${__dir}/convergence_cur_arr.py $1 $2 $3 $4 $5 $6 $7
#python ${__dir}/sCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7
#python ${__dir}/nCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7
#python ${__dir}/run_DMRG_pinning.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
python ./run_DMRG_pinning.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
#python ${__dir}/sCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
#python ./sCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
#python ${__dir}/nCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10
#python ./nCorr-func-GS.py $1 $2 $3 $4 $5 $6 $7 $8 $9 $10