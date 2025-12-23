#!/bin/bash
#SBATCH --job-name="download_betas_${SLURM_ARRAY_TASK_ID}"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-253%4
#SBATCH --mem=40G
#SBATCH --output=/jmsh/projects/researchers/bins14/AG_Israel/slurm_logs/200bp_bed_%A_%a.out
#SBATCH --error=/jmsh/projects/researchers/bins14/AG_Israel/slurm_logs/200bp_bed_%A_%a.err


# Specify the path to the config file
config="/jmsh/projects/researchers/bins14/AG_Israel/src/config.txt"

# Extract the sample name for the current $SLURM_ARRAY_TASK_ID and its variables 
ftp=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)

sample=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# output variables
wgbs_tools="/jmsh/projects/researchers/bins14/softwares/wgbs_tools/wgbstools"
bed_output="/jmsh/external/nihit/Israeli_methylation_dataset/bed"
beta_output="/jmsh/external/nihit/Israeli_methylation_dataset/beta"
output=${beta_output}/${sample}.hg38.beta
output_200bp="/jmsh/projects/researchers/bins14/AG_Israel/200bp_sum/"

<< 'COMMENT'
echo -e "Downloading beta file for ${sample}"



echo -e "wget -O ${output} ${ftp}"

wget -O ${output} ${ftp}
source /jmsh/projects/conda/minconda/bin/activate /jmsh/projects/conda/minconda/envs/wgbs_env



python3 ${wgbs_tools} view --genome hg38   ${beta_output}/${sample}.hg38.beta > ${bed_output}/${sample}.hg38.bed

COMMENT

counts="/jmsh/external/nihit/Israeli_methylation_dataset/bed/counts.bed"
function="mean"
function="sum"



echo -e "Generating 200bp file: ${output_200bp}/${sample}.200bp.hg38.bed"
source /jmsh/projects/conda/minconda/bin/activate /jmsh/projects/conda/minconda/envs/wgbs_env
bedtools  intersect -a ${counts} -b ${bed_output}/${sample}.hg38.bed -wao \
    |awk -vOFS="\t" '{print $1,$2,$3, ($6==-1?0:$7), ($6==-1?0:$8)}' \
    |bedtools merge -i - -d -1 -c 4,5 -o ${function} \
    |awk -vOFS="\t" '{if ($4>$5) $4=$5; print $1,$2,$3,sprintf("%.0f", $4),sprintf("%.0f", $5)}' \
    > ${output_200bp}/${sample}.${function}.200bp.hg38.bed