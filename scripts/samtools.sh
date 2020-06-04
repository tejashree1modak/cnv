#!/bin/sh
#BATCH --job-name="samtools"                                                                                                                                 ##SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="samtools_out.%A-%a"
#SBATCH --error="samtools_out.%A-%a"

set -e
echo "START" $(date)

module load bio/SAMtools/1.9-foss-2018b

infile=/data3/schwartzlab/tejashree/oyster_cnv/annot/gimap_genes.bed
outfile=/data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads.bam

echo "Infile ${infile}"
echo "Outfile ${outfile}"
echo "Starting samtools "


#samtools view -b -L "${infile}" /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.bam > /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads.bam

#samtools view -c /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads.bam #51054 reads mapped  

#samtools index -b /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads.bam /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads.bam.bai

#samtools flagstat "${outfile}"

#samtools view -c -f 0x8 "${outfile}" #This resulted in 156 which are the singletons per the flagstat output
#samtools view "${outfile}" > /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads_nohead.sam
#samtools view -b /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads_discordant.sam > /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads_discordant.bam

samtools index -b /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads_discordant.bam /data3/schwartzlab/tejashree/oyster_cnv/copy_loc/CL_1_gimap_reads_discordant.bam.bai

#echo "Outfile made ${outfile}"
echo "STOP" $(date)
echo "Done"
#(END)

