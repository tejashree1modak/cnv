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
#module load Python/3.6.4-foss-2018a
module load Pysam/0.15.1-foss-2018b-Python-2.7.15
module load numpy/1.11.0-foss-2017a-Python-2.7.13
module load svtyper/0.6.1-foss-2017a-Python-2.7.13
module load BCFtools/1.6-foss-2016b
module load SAMtools/1.5-foss-2017a
module load svtools/0.5.1-foss-2018b-Python-2.7.15
module load vawk/0.0.1-foss-2018b-Python-2.7.15

echo "modules loaded"

# Support for arguments passed in on the command line
# Run as sbatch --export input=filename .... 
#   if [ ! -f "$input" ];  then
#       echo "ERROR: '$input' not a valid file"
#       exit 1
#   fi
#   
#   case "$input" in
#   *.bam)  ;;  # all good
#   *)  echo "ERROR: '$input' has to end in *.bam"
#       exit 1
#       ;;
#   esac
#   
#   cd $(dirname $input)
#
#
cd /data3/schwartzlab/tejashree/bam/filtered/
#
#  echo "Starting samtools"
#  
#  samtools view -b -F 1294 "$input" > "${input%.bam}.discordants.unsorted.bam"
#  
#  samtools view -h "$input" \
#          | $EBROOTLUMPY/scripts/extractSplitReads_BwaMem -i stdin \
#              | samtools view -Sb - \
#                  > "${input%.bam}.splitters.unsorted.bam"
#  
#  # Sort both alignments
#  echo "Sorting with samtools"
#  file_prefix=$(basename "${input%.bam}")
#  samtools sort -T "${file_prefix}.dis" -o "${file_prefix}.discordants" "${file_prefix}.discordants.unsorted.bam"
#  samtools sort -T "${file_prefix}.split" -o "${file_prefix}.splitters" "${file_prefix}.splitters.unsorted.bam"
#
##Check sorting 
#  echo "Checking sorting"
#  python $EBROOTLUMPY/scripts/check_sorting.py \
#          "${file_prefix}.discordants.unsorted.bam" \
#              "${file_prefix}.splitters.unsorted.bam"
#  
#Run LUMPY (did not pass -P option,check if you need to)
#  echo "Running lumpy"
#     lumpyexpress \
#         -B "$input" \
#             -S "${input%.bam}.splitters" \
#                 -D "${input%.bam}.discordants" \
#                 -P \
#                     -o "${input%.bam}.vcf"
# 
###Run SVtyper 
#  echo "Running SVtyper"
#  svtyper \
#      -B "$input" \
#      -S "${input%.bam}.splitters" \
#      -i "${input%.bam}.vcf" \
#       > "${input%.bam}.gt.vcf" \

 # Sort and merge individual vcf files

   #python $EBROOTLUMPY/scripts/l_sort.py CL_1.F.vcf,CL_2.F.vcf,CL_3.F.vcf > CL_123_sorted.vcf

 # Sort vcf files 
#bcftools sort CL_1.F.gt.vcf -o CL_1.F.gt_sorted.vcf

 #zip vcf files
   #ls *.gt.vcf | xargs -n1 bgzip  
   #bgzip CL_1.F.gt.vcf 
   #bgzip CL_2.F.gt.vcf 
   #bgzip CL_3.F.gt.vcf 
 #index zipped files
 #bcftools index CL_1.F.gt.vcf.gz
    
 # Using bcftools to merge gt.vcf files from multiple samples 
   #bcftools merge -m id -O v -o 123merged.vcf *.gt.vcf.gz

 # Run svtools for next steps
# svtools lsort CL_1.F.gt.vcf CL_2.F.gt.vcf CL_3.F.gt.vcf \
# | bgzip -c > sorted.vcf.gz

#zcat sorted.vcf.gz \
#| svtools lmerge -i /dev/stdin -f 20 \
#| bgzip -c > merged.vcf.gz

mkdir -p gt
 zcat merged.vcf.gz \
| vawk --header '{  $6="."; print }' \
| svtools genotype \
  -B CL_1.F.bam \
  -l CL_1.F.bam.json \
| sed 's/PR...=[0-9\.e,-]*\(;\)\{0,1\}\(\t\)\{0,1\}/\2/g' - \
> gt/CL_1.F.vcf

echo "STOP" $(date)
echo "Done"



