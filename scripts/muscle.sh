#!/bin/bash
#SBATCH --job-name="muscle"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="muscle_out.%A-%a"
#SBATCH --error="muscle_out.%A-%a"

set -e
echo "START" $(date)

module load muscle/3.8.31

#muscle -in /data3/schwartzlab/tejashree/oyster_cnv/ifi44_dup_chr9.fa -out /data3/schwartzlab/tejashree/oyster_cnv/ifi44_dup_chr9.afa
muscle -in /data3/schwartzlab/tejashree/oyster_cnv/gimap_genes.fa -out /data3/schwartzlab/tejashree/oyster_cnv/gimap_genes.afa

echo "STOP" $(date)
echo "Done"
(END)
