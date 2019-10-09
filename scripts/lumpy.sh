#!/bin/sh
#SBATCH --job-name="lumpy"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="lumpy_out.%A-%a"
#SBATCH --error="lumpy_out.%A-%a"

set -e
echo "START" $(date)

module load LUMPY/0.3.0-foss-2018b
#module load LUMPY/0.2.13-foss-2016b
module load bio/SAMtools/1.9-foss-2018b
module load Python/2.7.15-foss-2018b
module load Pysam/0.15.1-foss-2018b-Python-2.7.15
module load numpy/1.11.0-foss-2017a-Python-2.7.13
module load svtyper/0.6.1-foss-2017a-Python-2.7.13

echo "modules loaded"

cd /data3/schwartzlab/tejashree/bam/filtered/

#echo "Starting samtools"

#samtools view -b -F 1294 /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.bam > /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.discordants.unsorted.bam

#samtools view -h /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.bam \
#        | $EBROOTLUMPY/scripts/extractSplitReads_BwaMem -i stdin \
#            | samtools view -Sb - \
#                > /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.splitters.unsorted.bam

# Sort both alignments
#echo "Sorting with samtools"
#samtools sort -T CL_1.F.dis -o CL_1.F.discordants CL_1.F.discordants.unsorted.bam
#samtools sort -T CL_1.F.split -o CL_1.F.splitters CL_1.F.splitters.unsorted.bam

#Check sorting 
#echo "Checking sorting"
#python $EBROOTLUMPY/scripts/check_sorting.py \
#        CL_1.F.discordants.unsorted.bam \
#            CL_1.F.splitters.unsorted.bam

#Run LUMPY
#echo "Running lumpy"
#    lumpyexpress \
#        -B /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.bam \
#            -S /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.splitters \
#                -D /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.discordants \
#                    -o /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.vcf

#Run SVtyper 
echo "Running SVtyper"
svtyper \
    -B /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.bam \
    -S data3/schwartzlab/tejashree/bam/filtered/CL_1.F.splitters \
    -i /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.vcf \
    > /data3/schwartzlab/tejashree/bam/filtered/CL_1.F.gt.vcf

echo "STOP" $(date)
echo "Done"



