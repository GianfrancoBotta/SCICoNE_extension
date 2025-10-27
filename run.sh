#!/bin/bash
#SBATCH --job-name=snakemake_qc
#SBATCH --ntasks=1      
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=30G     
#SBATCH --time=96:00:00      
#SBATCH --output snakemake_qc.log
#SBATCH --mail-type=END
#SBATCH --mail-user=gbotta@ethz.ch

source ~/.bashrc
conda activate snakemake
module load eth_proxy
snakemake --use-conda --use-singularity --cores 10 --singularity-args '-B /cluster/work/bewi/members/jgawron -B /scratch -B /cluster/work/bewi/members/gbotta:/cluster/work/bewi/members/gbotta:rw' --rerun-incomplete