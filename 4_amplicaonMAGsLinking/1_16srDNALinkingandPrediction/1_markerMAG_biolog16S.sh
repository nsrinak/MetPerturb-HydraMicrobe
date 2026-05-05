#!/bin/bash

#SBATCH --job-name=linking_filtered_biolog
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=500000
#SBATCH --time=48:00:00
#SBATCH --output=linking_filtered_biolog.out
#SBATCH --error=linking_filtered_biolog.err
#SBATCH --partition=base
#SBATCH --qos=normal

dir="."
dir_16s="./biolog_16s"

# Here we used both dereplicated bins and isolated genomes.

bins="./iso_drep_70_20"

# We combined sequencing reads from every samples (preprocessed files) into the same fasta file.
# So we have unmapped_all_metagenome_R1.fasta and unmapped_all_metagenome_R2.fasta as in command below.

MarkerMAG link -min_16s_len 200 -no_polish -no_cluster -p "unmapped_all_metagenome_filteted_biolog_16s" -r1 "${dir}/unmapped_all_metagenome_R1.fasta" -r2 "${dir}/unmapped_all_metagenome_R2.fasta" -marker "${dir_16s}/filtered_biolog_16s.fasta" -mag "${bins}" -x fa -t 32
