#!/bin/bash
#SBATCH --job-name="delly_Bam_files"
#SBATCH --time=999:00:00  # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=20   # processor core(s) per node
#SBATCH --mail-user="aakintomide@my.uri.edu"
#SBATCH --mail-type=END,FAIL
#SBATCH --output="out_oyster_delly"
#SBATCH --error="out_oyster_delly"

module load  delly/0.7.8-foss-2016b

cd $SLURM_SUBMIT_DIR
module load SAMtools/1.5-foss-2017a
module load BCFtools/1.6-foss-2016b

delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_1.bcf CL_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_2.bcf CL_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_3.bcf CL_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_4.bcf CL_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_5.bcf CL_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_6.bcf CL_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_1.bcf CLP_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_2.bcf CLP_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_3.bcf CLP_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_4.bcf CLP_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_5.bcf CLP_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CLP_6.bcf CLP_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_1.bcf CS_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_2.bcf CS_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_3.bcf CS_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_5.bcf CS_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_6.bcf CS_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CS_7.bcf CS_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_1.bcf DEBY_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_2.bcf DEBY_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_3.bcf DEBY_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_4.bcf DEBY_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_5.bcf DEBY_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o DEBY_6.bcf DEBY_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_1.bcf HC_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_3.bcf HC_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_4.bcf HC_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_5.bcf HC_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_6.bcf HC_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_7.bcf HC_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_1.bcf HC_VA_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_2.bcf HC_VA_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_3.bcf HC_VA_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_4.bcf HC_VA_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_5.bcf HC_VA_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HC_VA_6.bcf HC_VA_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HG_HG0F2.bcf HG_HG0F2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HG_HG2F1.bcf HG_HG2F1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HG_HG2M5.bcf HG_HG2M5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_1.bcf HI_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_2.bcf HI_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_3.bcf HI_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_4.bcf HI_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_5.bcf HI_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o HI_6.bcf HI_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LM_1_pool.bcf LM_1_pool.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LM_3.bcf LM_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LM_4.bcf LM_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LM_7.bcf LM_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LM_8.bcf LM_8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_1.bcf LOLA_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_2.bcf LOLA_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_3.bcf LOLA_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_4.bcf LOLA_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_5.bcf LOLA_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o LOLA_6.bcf LOLA_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_1.bcf NEH_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_2.bcf NEH_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_3.bcf NEH_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_4.bcf NEH_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_5.bcf NEH_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NEH_6.bcf NEH_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NG_NH0H4.bcf NG_NH0H4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NG_NH2F6.bcf NG_NH2F6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NG_NH2F8.bcf NG_NH2F8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o NG_NH2M1.bcf NG_NH2M1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_1.bcf OBOYS2_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_2.bcf OBOYS2_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_3.bcf OBOYS2_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_4.bcf OBOYS2_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_5.bcf OBOYS2_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o OBOYS2_6.bcf OBOYS2_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_1.bcf SL_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_2.bcf SL_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_3.bcf SL_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_4.bcf SL_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_5.bcf SL_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SL_6.bcf SL_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_10.bcf SM_10.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_11.bcf SM_11.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_12.bcf SM_12.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_7.bcf SM_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_8.bcf SM_8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o SM_9.bcf SM_9.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_1.bcf UMFS_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_2.bcf UMFS_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_3.bcf UMFS_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_4.bcf UMFS_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_5.bcf UMFS_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o UMFS_6.bcf UMFS_6.F.bam




delly merge -o sites.bcf *.bcf 


delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_1.geno.bcf CL_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_2.geno.bcf CL_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_3.geno.bcf CL_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_4.geno.bcf CL_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_5.geno.bcf CL_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CL_6.geno.bcf CL_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_1.geno.bcf CLP_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_2.geno.bcf CLP_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_3.geno.bcf CLP_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_4.geno.bcf CLP_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_5.geno.bcf CLP_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CLP_6.geno.bcf CLP_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_1.geno.bcf CS_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_2.geno.bcf CS_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_3.geno.bcf CS_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_5.geno.bcf CS_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_6.geno.bcf CS_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o CS_7.geno.bcf CS_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_1.geno.bcf DEBY_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_2.geno.bcf DEBY_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_3.geno.bcf DEBY_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_4.geno.bcf DEBY_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_5.geno.bcf DEBY_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o DEBY_6.geno.bcf DEBY_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_1.geno.bcf HC_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_3.geno.bcf HC_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_4.geno.bcf HC_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_5.geno.bcf HC_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_6.geno.bcf HC_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_7.geno.bcf HC_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_1.geno.bcf HC_VA_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_2.geno.bcf HC_VA_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_3.geno.bcf HC_VA_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_4.geno.bcf HC_VA_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_5.geno.bcf HC_VA_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HC_VA_6.geno.bcf HC_VA_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HG_HG0F2.geno.bcf HG_HG0F2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HG_HG2F1.geno.bcf HG_HG2F1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HG_HG2M5.geno.bcf HG_HG2M5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_1.geno.bcf HI_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_2.geno.bcf HI_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_3.geno.bcf HI_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_4.geno.bcf HI_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_5.geno.bcf HI_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o HI_6.geno.bcf HI_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LM_1_pool.geno.bcf LM_1_pool.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LM_3.geno.bcf LM_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LM_4.geno.bcf LM_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LM_7.geno.bcf LM_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LM_8.geno.bcf LM_8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_1.geno.bcf LOLA_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_2.geno.bcf LOLA_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_3.geno.bcf LOLA_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_4.geno.bcf LOLA_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_5.geno.bcf LOLA_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o LOLA_6.geno.bcf LOLA_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_1.geno.bcf NEH_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_2.geno.bcf NEH_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_3.geno.bcf NEH_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_4.geno.bcf NEH_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_5.geno.bcf NEH_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NEH_6.geno.bcf NEH_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NG_NH0H4.geno.bcf NG_NH0H4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NG_NH2F6.geno.bcf NG_NH2F6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NG_NH2F8.geno.bcf NG_NH2F8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o NG_NH2M1.geno.bcf NG_NH2M1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_1.geno.bcf OBOYS2_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_2.geno.bcf OBOYS2_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_3.geno.bcf OBOYS2_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_4.geno.bcf OBOYS2_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_5.geno.bcf OBOYS2_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o OBOYS2_6.geno.bcf OBOYS2_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_1.geno.bcf SL_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_2.geno.bcf SL_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_3.geno.bcf SL_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_4.geno.bcf SL_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_5.geno.bcf SL_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SL_6.geno.bcf SL_6.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_10.geno.bcf SM_10.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_11.geno.bcf SM_11.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_12.geno.bcf SM_12.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_7.geno.bcf SM_7.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_8.geno.bcf SM_8.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o SM_9.geno.bcf SM_9.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_1.geno.bcf UMFS_1.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_2.geno.bcf UMFS_2.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_3.geno.bcf UMFS_3.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_4.geno.bcf UMFS_4.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_5.geno.bcf UMFS_5.F.bam
delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -v sites.bcf -o UMFS_6.geno.bcf UMFS_6.F.bam

bcftools merge -m id -O b -o merged.bcf *.geno.bcf

delly filter -f germline -o germline.bcf merged.bcf


