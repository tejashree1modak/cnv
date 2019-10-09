#!/bin/bash
#SBATCH --job-name="delly_Bam_files"
#SBATCH --time=999:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="aakintomide@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="out_oyster_delly_2"
#SBATCH --error="out_oyster_delly_2"

module load  delly/0.7.8-foss-2016b

cd $SLURM_SUBMIT_DIR
module load SAMtools/1.5-foss-2017a
module load BCFtools/1.6-foss-2016b

bcftools index merged.bcf 
delly filter -f germline -o germline.bcf merged.bcf


