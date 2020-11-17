#!/bin/bash
#SBATCH --job-name="delly"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="delly_out.%A-%a"
#SBATCH --error="delly_out.%A-%a"

set -e
echo "START" $(date)

if [ -n "$V" ]; then
    # this is a dry run
    V=echo
fi

#Run script like this: 
# sbatch -a

#This script finds and maps SVs in samples by comparing sample genomes to reference using delly
#It takes inputs as reference genome fasta file and BAM files of samples 
#To run script:

module load  delly/0.7.8-foss-2016b
module load SAMtools/1.5-foss-2017a
module load BCFtools/1.6-foss-2016b
module load VCFtools/0.1.16-foss-2019a-Perl-5.28.1
module load tabix/0.2.6-GCCcore-7.3.0

echo "modules loaded"

reference=/data3/schwartzlab/tejashree/masked_genome/masked.cvir.genome.fasta
delly_out=/data3/schwartzlab/tejashree/masked_genome/delly_out

cd $delly_out

FILES=(
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CL_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CLP_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/CS_7.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/DEBY_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_7.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HC_VA_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HG_HG0F2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HG_HG2F1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HG_HG2M5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/HI_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LM_1_pool.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LM_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LM_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LM_7.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LM_8.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/LOLA_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NEH_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NG_NH0H4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NG_NH2F6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NG_NH2F8.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/NG_NH2M1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/OBOYS2_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SL_6.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_10.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_11.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_12.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_7.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_8.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/SM_9.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_1.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_2.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_3.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_4.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_5.F.bam
/data3/schwartzlab/tejashree/masked_genome/masked/filtered/UMFS_6.F.bam
)

if [ -n "$SLURM_ARRAY_TASK_ID" ] ; then
    # this is a an array job
    # SLURM_ARRAY_TASK_ID will determine which file among this will is to be processed
    if [ "$SLURM_ARRAY_TASK_ID" -gt ${#FILES[*]} ]; then
        echo "Index $SLURM_ARRAY_TASK_ID more than ${#FILES[*]}"
        exit 1
    fi

    let SLURM_ARRAY_TASK_ID-- # SLURM_ARRAY_TASK_ID is 1 based, arrays are 0 based
    infile=${FILES[$SLURM_ARRAY_TASK_ID]}

    if [ ! -f "$infile" ]; then
        echo "$infile not found, aborting"
        exit 1
    fi
fi 

if [ "$STEP" -eq 1 ]; then
    # check if index file is there, else create an index
    if [ ! -f "${infile}.bai" ] ; then
        echo "Creating index for ${infile} (`date`)"
        $V samtools index "$infile"
        echo "Done index for ${infile} (`date`)"
    fi

    out_bfc=$(basename $infile | sed -e 's/.F.bam//g')
    $V delly call -g "${reference}" -o "${out_bfc}.bfc"  "$infile"
elif [ "$STEP" -eq 2 ]; then
    $V delly merge -o sites.bcf *.bfc
elif [ "$STEP" -eq 3 ]; then
    out_bfc=$(basename $infile | sed -e 's/.F.bam//g')
    sites="sites.bcf"

    if [ ! -f "$sites" ]; then
        echo "$sites file not found in $PWD"
        exit 1
    fi
   $V delly call -g "${reference}" -v "$sites" -o "${out_bfc}.geno.bcf" "$infile"
elif [ "$STEP" -eq 4 ]; then
    out_file=merged.bcf
    $V bcftools merge -m id -O b -o $out_file *.geno.bcf
elif [ "$STEP" -eq 5 ] ; then
    out_file=merged.bcf
    $V bcftools index $out_file
    $V delly filter -f germline -o germline.bcf $out_file
elif [ "$STEP" -eq 6 ] ; then
    $V bcftools view -O v germline.bcf > germline.vcf 
    #Convert this vcf to tsv by opening in vim and delete the header saved as germline.nohead.vcf 
elif [ "$STEP" -eq 7 ] ; then
    #$V bgzip germline.vcf
    #$V tabix -p vcf germline.vcf.gz
    #$V zcat germline.vcf.gz | vcf-to-tab > germline.tab
else
    echo "STEP variable not set"
    exit 1
fi

echo "STOP" $(date)
