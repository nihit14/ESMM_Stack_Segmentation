#### Prepare the chromsome wise matrix
For each chromosome from bed file which is converted from beta file using wgbstools 
Here is the discription of beta file 
beta file is the simplest representation of methylation data. It is a binary file with a fixed size of 2 * 8 * NR_SITES bytes (~54MB for hg19), holding a matrix of uint8 values with dimensions of (NR_SITES x 2). For each of the NR_SITES CpG sites in the genome, it holds 2 values: the #meth and #covered. Meaning, the i'th row in the matrix corresponds to the i'th CpG site:

meth: the number of times the i'th site (CpGi) site is observed in a methylated state.
coverage: the total number of times i'th site (CpGi) is observed. #coverage==0 is equivalent to a missing value (NaN).
CpGi beta value is obtained by dividing #meth / #coverage.

# Read all the files and then use divide for each chromosome
cd /jmsh/external/nihit/Israeli_methylation_dataset/matrices 

# create folder structure 
for chromosomes in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
do 
    mkdir chr${chromosomes}
done

# the following command will go through each file and divide based on column1 which is chromosome and for each file it extracts the file name and save it in the folder ${chromosome}/${samplename}_{chromosome}.bed 
# only first three columns are saved # change to print $5 incase you need coverage and $4 if you need methylation (Keep in mind here 5th column is coverage and 4th column is methyated reads only because we get it from beta files)
# NR > 1 in case you have file with header for below command. 

for files in $(ls ../bed/*bed |grep -v counts )
do 
    fname=$(basename $files .hg38.bed)
    awk -v a="$fname" -vOFS="\t" '{print $1,$2,$3 > $1"/"a"_"$1".bed"}' $files 
done

# above command was cancelled as I could think of efficient way to get all coordinates across files I simply need to extract chromosome wise second column from each file and this can be used as the coordinate # Adapted command only print $2
for files in $(ls ../bed/*bed |grep -v counts ); do fname=$(basename $files .hg38.bed);awk -v a="$fname" -vOFS="\t" '{print $2 > $1"/"a"_"$1".txt"}' $files ;done

# Then cat all txt files chromosome wise |sort them and find uniq entries.  
for chromosomes in $(ls c* -d ) ;do echo $chromosomes;cat ${chromosomes}/*txt |sort |uniq > ${chromosomes}.txt  ;done


# create proper bedfile for the all CpG coordinates
for files in chr*txt ; do fname=$(basename $files .txt); sort -k1,1n $files |awk -v a="$fname" -vOFS="\t" '{print a,$1,$1+2 }' > ${fname}_sorted.bed  ;done


# smoothing of all samples for 200 and 1000 bp 
ls /jmsh/external/nihit/Israeli_methylation_dataset/bed/*.hg38.bed |nl |awk -vOFS="\t" 'BEGIN{print "ArrayTaskID\tSample\tSmoothing_window"}{print $0,200}' > config_smoothing.txt

ls /jmsh/external/nihit/Israeli_methylation_dataset/bed/*.hg38.bed |awk -vOFS="\t" 'BEGIN{i=253}{i+=1;print i,$0,1000}' >> config_smoothing.txt

# First only run for 200bp




# Since R wasn't able to handle the large files we used bedtools to merge all files to get column wise Cov and Meth for each file. 
# working directory 
cd /jmsh/external/nihit/Israeli_methylation_dataset/matrices_smooth_200bp

# make chromosome wise bed files 
for files in $(ls ../bed_smooth_200bp/*_200_smooth.txt |grep -v counts )
do 
    fname=$(basename $files _200_smooth.txt)
    echo $fname 
    awk -v a="$fname" -vOFS="\t" 'NR>1{print $0 > "chr"$1"/"a"_chr"$1".bed"}' $files
done

# alternative to above command. Using find and one can paralllize it using xargs 
find  ../bed_smooth_1000bp/ -name "*_1000_smooth.txt"  \
    |  xargs -n 1 sh -c 'files="$1"; fname=$(basename $files _1000_smooth.txt) ; echo ${fname}; awk -v a="${fname}" -vOFS="\t" '\''NR>1{print $0 > "chr"$1"/"a"_chr"$1".bed"}'\'' ${files}' _

# first list all 253 files and keep this file as this is the order that will be also used later on for yaml file as well as for annotation. 
ls chr1/*bed  > files_bed.txt

# create initial bed file with only coordinates from that particular chromosome
cp chr1_sorted.bed test.bed

# command that will iterate over all the files and create the final bed file which will be test.bed

# for this first create intially all the files with coordinates 
for chromosomes in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M ; do cp chr${chromosomes}_sorted.bed chr${chromosomes}_sorted_final.txt ;done

# improve the same loop now to use it for the creating and replacing the file after each intersect 
# I keep extra the files_bed.txt to keep the order consistent 

for chromosomes in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y M
do 
    while read bedfiles
    do
        fname=$(basename $bedfiles _chr1.temp.bed)
        bed_to_intersect="chr${chromosomes}/${fname}_chr${chromosomes}.temp.bed"
        echo -e "Starting chromosome ${chromosomes} and sample ${fname}. " 
        bedtools intersect -a chr${chromosomes}_sorted_final.txt -b  $bed_to_intersect -wao \
        |awk -vOFS="\t" '{$(NF-5)=$(NF-4)=$(NF-3)=$(NF)=""; print $0}' \
        |sed 's/\t\+/\t/g;s/^\t//' \
        |sed 's/\t$//' \
        > temp.bed && mv temp.bed chr${chromosomes}_sorted_final.txt
    done < files_bed.txt
    echo -e "Finished chromosome ${chromosomes}."     
done > loop.log 2>&1


# Now we remove (chr) for the raw matrix from all chromosomes and store all the matrix in organised folder where it is only read rights with metadata
# create header for file 
header=$(for seq in $(seq 1 253); do echo -e 'Cov\tMeth' ;done |tr '\n' '\t')

# mkdir 
mkdir /jmsh/external/nihit/Israeli_methylation_dataset/matrices/{raw,smooth_200bp,smooth_1000bp}

# remove chr for raw matrix and store them in right folder 
cd /jmsh/external/nihit/Israeli_methylation_dataset/matrices
for files in $(ls *sorted_final.txt)
do 
    fname=$(basename $files _sorted_final.txt)
    echo $fname
    awk -vOFS="\t" -v a="${header}" 'BEGIN{print "chr\tstart\tend",a}{gsub(/chr/,"");print $0}' $files \
    |sed 's/\t$//'  \
    > raw/${fname}.txt
done


# creating methlyation ratio matrices 
cd /jmsh/external/nihit/Israeli_methylation_dataset/matrices
# for each chromosome find the corresponding chr file and make a methylation ratio 
for files in $(ls */chr*txt |grep -v old) ; do echo -e "Processing ${files}" ;awk -vOFS="\t" '{printf "%s%s%s", $1,OFS,$2; printf OFS $3; for(i=4;i<NF;i+=2) printf OFS int(($(i)==0)?0:($(i+1)/$i)*100); printf ORS}' $files > ../methylation_matrices/${files} ;done

# concatenate each file to get complete matrix
cat  ../methylation_matrices/smooth_1000bp/chr*txt |grep -v nan > ../methylation_matrices/temp.txt && mv ../methylation_matrices/temp.txt ../methylation_matrices/raw_methylation_matrix.txt

# extract header and add the header to the main file 
header=$(cat ../methylation_matrices/column_header.txt |cut -f 1 -d '_' |tr '\n' '\t')
awk -vOFS="\t" -v a="$header" 'BEGIN{sub(/\t$/, "", a); print a}{print $0}' ../methylation_matrices/methylation_matrix_smooth-1000bp.txt  > ../methylation_matrices/temp.txt && mv ../methylation_matrices/temp.txt ../methylation_matrices/methylation_matrix_smooth-1000bp.txt