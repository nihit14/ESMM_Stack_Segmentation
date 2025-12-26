# default values 
binsize=200
shift_bp=75
paired_end=ignore
cores=100
smoothing_window=1000 # ? change according to research questions (LMR ~ 200bp (n=25) PMDs 1000bp (n=75))

while getopts ":hi:f:o:b:" OPTION
do
    case $OPTION in
                i) 
                    input=$OPTARG ;; 
                f)
                    bed_file=$OPTARG ;;
                o)
                    output=$OPTARG ;;
                b) 
                    binsize=$OPTARG ;;
    esac
done

smoothing_window=${binsize}
filename="${output}/${input}"
methyl_bed=${bed_file}


smoothing_script_R="/jmsh/projects/researchers/bins14/AG_Israel/src/smoothing/smoothing.r"

echo -e "Rscript --vanilla ${smoothing_script_R} ${methyl_bed} ${smoothing_window} ${filename} ${cores}"

Rscript --vanilla ${smoothing_script_R} ${methyl_bed} ${smoothing_window} ${filename} ${cores}
