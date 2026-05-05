#!/bin/bash

#SBATCH --job-name=linking_metagenomic_matam
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500000
#SBATCH --time=48:00:00
#SBATCH --output=linking_metagenomic_matam.out
#SBATCH --error=linking_metagenomic_matam.err
#SBATCH --partition=base
#SBATCH --qos=normal

dir="."

# I moved all matam 16S rRNA results (99% similarity results) to the same directory for convenience

dir_16s="./matam_16s/16s_99_similarity"
bins="./iso_drep_70_20"

MarkerMAG link -p "unmapped_all_metagenome" -r1 "${dir}/unmapped_all_metagenome_R1.fasta" -r2 "${dir}/unmapped_all_metagenome_R2.fasta" -marker "${dir_16s}/test_110_samples.fasta" -mag "${bins}" -x fa -t 32

