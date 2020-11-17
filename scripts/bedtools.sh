#!/bin/sh
#SBATCH --job-name="bedtools"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="bedtools_out.%A-%a"
#SBATCH --error="bedtools_out.%A-%a"

set -e
echo "START" $(date)

module load bio/BEDTools/2.26.0-foss-2016b 

#infile=/data3/schwartzlab/tejashree/oyster_cnv/gimap_genes.bed
#outfile=/data3/schwartzlab/tejashree/oyster_cnv/gimap_genes.fa

#echo "Infile ${infile}"
echo "Starting bedtools intersect"

#bedtools getfasta -fi /data3/marine_diseases_lab/tejashree/Bio_project_SRA/cvir_genome/cvir_edited.fa -bed ${infile} -fo ${outfile} -name

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oysterduplicate_sort.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Cvir_genome_repeats_mod.bed -wo

#bedtools merge -i /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Cvir_genome_repeats_mod_strand.bed -c 1,5 -o count,collapse > /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Cvir_genome_repeats_mod_strand_merged.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oysterduplicate_sort.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Cvir_repeats_merged.bed -wo > /data3/schwartzlab/tejashree/oyster_cnv/repeat/dup_repeat_merged_overlap.bed 

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/virginica_exp_Cvir_XP_BED_info_unique_shortened5.bed.txt -wo > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_cv_expanded_unique_overlap.bed 

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/genome_feat.bed -wo > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_genome_feat_overlap.bed 

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/exon_annot.bed -wo -f 1 > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_exon_full_overlap.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/gene_annot.bed -wo -f 1 > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_gene_full_overlap.bed

#bedtools merge -i /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -c 1,2,3 -o count,collapse,collapse  > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups_merged.bed

#bedtools merge -i /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups_LM.bed -c 1,2,3 -o count,collapse,collapse  > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups_LM_merged.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/intron_intergenic.bed -wo -f 1 > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_intron_intergenic_full_overlap.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intron_Merged.bed -wo -f 1 > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_intron_full_overlap.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intron_Merged.bed -wo > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_intron_overlap.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intergenic_Merged.bed -wo -f 1 > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_intergenic_full_overlap.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intergenic_Merged.bed -wo  > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_intergenic_overlap.bed

## intersect between all dups and haplotigs
#bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/haplotigs.bed -wo  > /data3/schwartzlab/tejashree/oyster_cnv/dups/filtration/dup_haplotigs_overlap.bed

## No selected lines dups merge
#bedtools merge -i /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups.bed -c 1,2,3 -o count,collapse,collapse  > /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups_merged.bed

## Intersect between No selected lines dups and genomic features
# Dups commpletely in EXONS 
bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/exon_annot.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/dup_nosel_exon_full_overlap.bed
# Dups commpletely in GENES
bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/gene_annot.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/dup_nosel_gene_full_overlap.bed
# Dups commpletely in INTRONS 
bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intron_Merged.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/dup_nosel_intron_full_overlap.bed
# Dups commpletely in INTERGENIC 
bedtools intersect -a /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/cvir_filtered_nosel_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intergenic_Merged.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/oyster_cnv/oyster_cnv_nosel/dup_nosel_intergenic_full_overlap.bed


#echo "Outfile made ${outfile}"
echo "STOP" $(date)
echo "Done"
#(END)
