#!/bin/bash

#SBATCH --job-name=muscle_1
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500000
#SBATCH --time=48:00:00
#SBATCH --output=muscle_1.out
#SBATCH --error=muscle_1.err
#SBATCH --partition=base
#SBATCH --qos=normal

# combined_all_16s.fasta contains 16S rRNA sequence from biolog experiment sequencing, matam prediction, barrnap prediction, and GTDBtk database

muscle -align combined_all_16s.fasta -output combined_all_16s_aligned.fasta
Gblocks combined_all_16s_aligned.fasta -t=d
FastTree -gtr -nt all_16s_alignment.fasta-gb > tree_16SBiolog_16SMatam_16SBarrnap_16SGTDB

