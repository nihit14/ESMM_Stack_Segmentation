#!/bin/bash
#SBATCH --job-name="smoothing_${SLURM_ARRAY_TASK_ID}"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --array=1-506%4
#SBATCH --mem=500G
#SBATCH --output=/jmsh/projects/researchers/bins14/slurm_logs/smoothing/smoothing_%A_%a.out
#SBATCH --error=/jmsh/projects/researchers/bins14/slurm_logs/smoothing/smoothing_%A_%a.err

# Specify the path to the config file
config="/jmsh/projects/researchers/bins14/AG_Israel/src/config_smoothing.txt"


# Extract the sample file and smoothing window for the current $SLURM_ARRAY_TASK_ID and its variables 
sample_bed=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
smoothing_window=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# extract the sample name from file
sample_name=$(basename $sample_bed .hg38.bed)
echo ${sample_name}

# set variables for scripts and output
counts_scripts="/jmsh/projects/researchers/bins14/AG_Israel/src/smoothing/wgbs_counts.sh"
smooth_container="/jmsh/projects/researchers/bins14/softwares/bssmooth/bsseq.sif"
output="/jmsh/external/nihit/Israeli_methylation_dataset/bed_smooth_${smoothing_window}bp"

echo "Running the smoothing container !!!"
echo -e "singularity exec -B /jmsh -B /walter ${smooth_container} bash ${counts_scripts} -i ${sample_name} -f ${sample_bed} -o ${output} -b ${smoothing_window}"

singularity exec -B /jmsh -B /walter ${smooth_container} bash ${counts_scripts} -i ${sample_name} -f ${sample_bed} -o ${output} -b ${smoothing_window}
