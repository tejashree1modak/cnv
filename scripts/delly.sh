BATCH --job-name="stringtie"
#SBATCH --time=9999:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --output="delly_out.%A-%a"
#SBATCH --error="delly_out.%A-%a"

set -e
echo "START" $(date)

#This script finds and maps SVs in samples by comparing sample genomes to reference using delly
#It takes inputs as reference genome fasta file and BAM files of samples 
#To run script:

module load  delly/0.7.8-foss-2016b
module load SAMtools/1.5-foss-2017a
module load BCFtools/1.6-foss-2016b
echo "modules loaded"



delly call -g /data3/schwartzlab/ayomide/Oyster_CNV/Genome/reference.fasta -o CL_1.bcf CL_1.F.bam
