#!/bin/bash
#SBATCH --array=1-2
#SBATCH --mem=200G
#SBATCH --cpus-per-task=18
module load TensorFlow/2.6.2-foss-2021a
source ~/virtualenvs/Dashain/bin/activate

BASEDIR=/home/t326h379/Test_Protein_Proteomics_123

export FILENAME=$(ls ${BASEDIR}/*.fasta | sed -n ${SLURM_ARRAY_TASK_ID}p)
python ~/Prot_T5_file_generation_of_Proteins.py
