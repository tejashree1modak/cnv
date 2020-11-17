#!/bin/sh
#SBATCH --job-name="masked_bedtools"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="masked_bedtools_out.%A-%a"
#SBATCH --error="masked_bedtools_out.%A-%a"

set -e
echo "START" $(date)

module load bio/BEDTools/2.26.0-foss-2016b

#infile=/data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup.bed

#bedtools sort -i ${infile} > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_sort.bed

#bedtools intersect -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_sort.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Cvir_repeats_merged.bed -wo > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_repeat_merged_overlap.bed

#bedtools merge -i /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_cvir_filtered_dups.bed  -c 1,2,3 -o count,collapse,collapse > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_merged_cvir_filtered_dups.bed

#bedtools intersect -wa -wb -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_sort.bed -b /data3/schwartzlab/tejashree/oyster_cnv/annot/Oyster_gene_sort.bed > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_gene_overlap.bed


## Intersect between dups obtained from MASKED genome  and genomic features
# Dups commpletely in EXONS
bedtools intersect -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/exon_annot.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_exon_full_overlap.bed
# Dups commpletely in GENES
bedtools intersect -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/gene_annot.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_gene_full_overlap.bed
# Dups commpletely in INTRONS
bedtools intersect -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intron_Merged.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_intron_full_overlap.bed
# Dups commpletely in INTERGENIC
bedtools intersect -a /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_cvir_filtered_dups.bed -b /data3/schwartzlab/tejashree/oyster_cnv/cvir_genome/Oyster_intergenic_Merged.bed -wo -f 1.0 > /data3/schwartzlab/tejashree/masked_genome/delly_out/masked_dup_intergenic_full_overlap.bed



echo "END" $(date)
END
