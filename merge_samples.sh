#!/bin/bash
#SBATCH --job-name="smoothing_${SLURM_ARRAY_TASK_ID}"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --array=1-25%4
#SBATCH --mem=500G
#SBATCH --output=/jmsh/projects/researchers/bins14/slurm_logs/merge_samples/merge_samples_%A_%a.out
#SBATCH --error=/jmsh/projects/researchers/bins14/slurm_logs/merge_samples/merge_samples_%A_%a.err

# Specify the path to the config file
config="/jmsh/projects/researchers/bins14/AG_Israel/src/merge_samples_config.txt"


# Extract the sample file and smoothing window for the current $SLURM_ARRAY_TASK_ID and its variables 
chromosomes=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
smoothing_window=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

echo -e "Smoothing window of ${smoothing_window} for chromosome ${chromosomes}"


# set the variables
chr_dir="/jmsh/external/nihit/Israeli_methylation_dataset/matrices/" 
working_dir="/jmsh/external/nihit/Israeli_methylation_dataset/matrices_smooth_${smoothing_window}bp"

# copy the chromosome coordinates as the starting point
txt_merged_file="${working_dir}/chr${chromosomes}_sorted_final.txt"
awk -vOFS="\t" '{gsub(/chr/,""); print }' ${chr_dir}chr${chromosomes}_sorted.bed > ${txt_merged_file}

# bedfile to read the order of the files and add one bedfile at a time 
bed_file_list="/jmsh/projects/researchers/bins14/AG_Israel/src/files_bed.txt"

# activate the conda environment 
source /jmsh/projects/conda/minconda/bin/activate /jmsh/projects/conda/minconda/envs/bins14_core/

# run the loop to add sequentially a bed file at a time
while read bedfiles
do
    fname=$(basename $bedfiles _chr1.temp.bed)
    bed_to_intersect="${working_dir}/chr${chromosomes}/${fname}_chr${chromosomes}.bed"
    echo -e "Intersect between ${txt_merged_file}, ${bed_to_intersect}" 
    bedtools intersect -a ${txt_merged_file} -b  $bed_to_intersect -wao \
    |awk -vOFS="\t" '{$(NF-5)=$(NF-4)=$(NF-3)=$(NF)=""; print $0}' \
    |sed 's/\t\+/\t/g;s/^\t//' \
    |sed 's/\t$//' \
    > ${working_dir}/chr${chromosomes}_temp.bed && mv ${working_dir}/chr${chromosomes}_temp.bed ${txt_merged_file}
done < ${bed_file_list} > ${working_dir}/chr${chromosomes}_loop.log 2>&1

echo -e "Finished chromosome ${chromosomes}."   