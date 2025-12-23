#### Prepare the chromsome wise matrix
For each chromosome from bed file which is converted from beta file using wgbstools 
Here is the discription of beta file 
beta file is the simplest representation of methylation data. It is a binary file with a fixed size of 2 * 8 * NR_SITES bytes (~54MB for hg19), holding a matrix of uint8 values with dimensions of (NR_SITES x 2). For each of the NR_SITES CpG sites in the genome, it holds 2 values: the #meth and #covered. Meaning, the i'th row in the matrix corresponds to the i'th CpG site:

meth: the number of times the i'th site (CpGi) site is observed in a methylated state.
coverage: the total number of times i'th site (CpGi) is observed. #coverage==0 is equivalent to a missing value (NaN).
CpGi beta value is obtained by dividing #meth / #coverage.

# Read all the files and then use divide for each chromosome
cd /jmsh/external/nihit/Israeli_methylation_dataset/matrices 

# the following command will go through each file and divide based on column1 which is chromosome and for each file it extracts the file name and save it in the folder ${chromosome}/${samplename}_{chromosome}.bed 
# only first three columns are saved # change to print $5 incase you need coverage and $4 if you need methylation (Keep in mind here 5th column is coverage and 4th column is methyated reads only because we get it from beta files)
for files in $(ls ../bed/*bed |grep -v counts ); do fname=$(basename $files .hg38.bed);awk -v a="$fname" -vOFS="\t" '{print $1,$2,$3 > $1"/"a"_"$1".bed"}' $files ;done

# above command was cancelled as I could think of efficient way to get all coordinates across files I simply need to extract chromosome wise second column from each file and this can be used as the coordinate # Adapted command only print $2
for files in $(ls ../bed/*bed |grep -v counts ); do fname=$(basename $files .hg38.bed);awk -v a="$fname" -vOFS="\t" '{print $2 > $1"/"a"_"$1".txt"}' $files ;done

# Then cat all txt files chromosome wise |sort them and find uniq entries.  
for chromosomes in $(ls c* -d ) ;do echo $chromosomes;cat ${chromosomes}/*txt |sort |uniq > ${chromosomes}.txt  ;done


# create proper bedfile for the all CpG coordinates
for files in chr*txt ; do fname=$(basename $files .txt); sort -k1,1n $files |awk -v a="$fname" -vOFS="\t" '{print a,$1,$1+2 }' > ${fname}_sorted.bed  ;done
