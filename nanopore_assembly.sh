#!/bin/bash



# Process and assemble Oxford nanopore sequencer data in fast5 format
# Include file format conversion, stats, qc and de novo assembly



#script version
version="0.4.7"


######################
#                    #
#    User Defined    #
#                    #
######################


# Analysis folder
export baseDir=/media/2TB_NVMe/NPWGS_20200217

# Reads
export fast5=/media/30tb_raid10/data/Mbovis/canada/nanopore/NPWGS_20200217/MBWGS183/20200218_2125_MN28835_ACA114_b6146520/fast5


# MinION sequencing kits
kit="SQK-LSK109"  # Library kit: SQK-LSK109, SQK-LSK108 (discontinued), SQK-RBK004 (rapid),
flow="FLO-MIN106"  # Flowcell kit: FLO-MIN106, FLO-FLG001
bk="EXP-NBD104"  # Barcoding kit: SQK-RBK004 (rapid), EXP-NBD104 (1-12), EXP-NBD114 (13-24)

# Barcode to sample name file
bc_desc=/media/30tb_raid10/data/meta90_low_parasites/minion/bc.txt

# Program location
export prog=""${HOME}"/prog"
export scripts=""${HOME}"/scripts"

# Maximum number of cores used per sample for parallel processing
# A highier value reduces the memory footprint.
export maxProc=4

# Estimated genome size in bp
# export size=4850000  # salmonella
# export size=6000000  # Pseudomonas
export size=4400000  # Mbovis

# Set assembler to use: unicycler, unicycler_hybrid, canu shasta, or flye
export assembler="shasta"

# Set smallest contig size for assemblies
export smallest_contig=1000

# #Kraken DB to use
# export kraken2_db="/media/30tb_raid10/db/kraken2/standard"  # refseq
export kraken2_db="/media/30tb_raid10/db/kraken2/nt"

#Annotation
export kingdom="Bacteria"
# export genus="Salmonella"
# export species="enterica"
export genus="Pseudomonas"
export species="syringaea"
export gram="neg"
export locus_tag="TOCHANGE"
export centre="OLF"


######################
#                    #
#     Resources      #
#                    #
######################


# Computer performance
export cpu=$(nproc) #total number of cores
export mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


#######################
#                     #
#   Data Stucture     #
#                     #
#######################


# Folder structure
export logs=""${baseDir}"/logs"
export qc=""${baseDir}"/qc"
export fastq=""${baseDir}"/fastq"
export trimmed=""${baseDir}"/trimmed"
export filtered=""${baseDir}"/filtered"
export basecalled=""${baseDir}"/basecalled"
export assemblies=""${baseDir}"/assemblies"
export polished=""${baseDir}"/polished"
export annotation=""${baseDir}"/annotation"
export phaster=""${baseDir}"/phaster"
export amr=""${baseDir}"/amr"


# Create folders if do not exist
# "||" if test is false
# "&&" if test is true
[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$qc" ] || mkdir -p "$qc"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$filtered" ] || mkdir -p "$filtered"
[ -d "$basecalled" ] || mkdir -p "$basecalled"
[ -d "$assemblies" ] || mkdir -p "$assemblies"
[ -d "$polished" ] || mkdir -p "$polished"
[ -d "$annotation" ] || mkdir -p "$annotation"
[ -d "$phaster" ] || mkdir -p "$phaster"
[ -d "$amr" ] || mkdir -p "$amr"


################
#              #
#     Log      #
#              #
################


# Date
echo -e "$(date)\n" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt

#script version
echo -e "\nnanopore_assembly.sh version "$version"\n" | tee -a "${logs}"/log.txt  # $0


####################
#                  #
#   Dependencies   #
#                  #
####################


# Check for dependencies and log versions

# Java
if hash java 2>/dev/null; then 
    java -version 2>&1 1>/dev/null | grep "java version" | tr -d '"' | tee -a "${logs}"/log.txt
else
    echo >&2 "java was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Python
if hash python3 2>/dev/null; then 
    python3 --version | tee -a "${logs}"/log.txt
else
    echo >&2 "python3 was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Guppy
if hash guppy_basecaller 2>/dev/null; then
    v=$(guppy_basecaller --version | rev | cut -d " " -f 1 | rev)
    echo "Guppy v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "Guppy was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# FastQC
if hash fastqc 2>/dev/null; then 
    fastqc -v | tee -a "${logs}"/log.txt
else
    echo >&2 "fastQC was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# MultiQC
source activate multiqc
if hash multiqc 2>/dev/null; then
    multiqc --version | tee -a "${logs}"/log.txt
    source deactivate
else
    echo >&2 "multiqc was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Kraken2
if hash kraken2 2>/dev/null; then
    kraken2 --version | grep -F "version" | tee -a "${logs}"/log.txt
else
    echo >&2 "kraken2 was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Porechop
if hash porechop 2>/dev/null; then
    v=$(porechop --version)
    echo "Porechop v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "porechop was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Filtlong
if hash filtlong 2>/dev/null; then
    filtlong --version | tee -a "${logs}"/log.txt
else
    echo >&2 "filtlong was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Flye
if hash flye 2>/dev/null; then
    v=$(flye --version)
    echo "flye v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "flye was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Unicycler
if hash unicycler 2>/dev/null; then
    unicycler --version | tee -a "${logs}"/log.txt
else
    echo >&2 "unicycler was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Shasta
if hash shasta 2>/dev/null; then
    shasta --version | grep -F "Shasta" | tee -a "${logs}"/log.txt
else
    echo >&2 "shasta was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Medaka
source activate medaka
if [ $? -eq 0 ]; then  # environment was found and activated without error
    medaka --version | tee -a "${logs}"/log.txt
    source deactivate
else
    echo >&2 "medaka was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

#QUAST
source activate quast
if hash quast 2>/dev/null; then
    quast --version | tee -a "${logs}"/log.txt
    source deactivate
else
    echo >&2 "QUAST was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Mummer
if hash porechop 2>/dev/null; then
    v=$(mummer --version)
    echo "Mummer v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "mummer was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Prokka
source activate prokka
if [ $? -eq 0 ]; then  # environment was found and activated without error
    prokka --version | tee -a "${logs}"/log.txt
    source deactivate
else
    echo >&2 "prokka was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# ResFinder
perl "${prog}"/resfinder/resfinder.pl -h
if [ $# -eq 0 ]; then
    version=$(perl "${prog}"/resfinder/resfinder.pl -h | grep "Current" \
        | cut -d ":" -f 2 | tr -d " ")
    echo "resfinder v"${version}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "resfinder was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# minimap2
if hash minimap2 2>/dev/null; then
    v=$(minimap2 --version 2>&1)
    echo "minimap2 "$v"" | tee -a "${logs}"/log.txt
else
    echo >&2 "minimap2 was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Samtools
if hash samtools 2>/dev/null; then
    samtools --version | grep -F 'samtools' | tee -a "${logs}"/log.txt
else
    echo >&2 "samtools was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Qualimap
if hash qualimap 2>/dev/null; then
    qualimap --version | grep -F "QualiMap" | tee -a "${logs}"/log.txt
else
    echo >&2 "qualimap was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# Blast
if hash blastn 2>/dev/null; then
    v=$(blastn -version | grep -F "blastn" | cut -d " " -f 2)
    echo "blastn v"${v}"" | tee -a "${logs}"/log.txt
else
    echo >&2 "blastn was not found. Aborting." | tee -a "${logs}"/log.txt
    exit 1
fi

# nanoQC
# Phaster
# CARD + RGI


####################################
#                                  #
#   Basecalling + Demultiplexing   #
#                                  #
####################################


guppy_bin_folder=$(dirname $(which guppy_basecaller))
guppy_data_folder="${guppy_bin_folder%/bin}/data"

#####   GPU mode   #####
#runners * chunks_per_runner * chunk_size < 100000 * [max GPU memory in GB]
# 2 * 1000 * 1000 < 100000 * 24

guppy_basecaller \
    --input_path "$fast5" \
    --save_path "$basecalled" \
    --kit "$kit" \
    --flowcell "$flow" \
    --recursive \
    --records_per_fastq 0 \
    --compress_fastq \
    --disable_pings \
    --model_file "${guppy_data_folder}"/template_r9.4.1_450bps_hac.jsn \
    --calib_detect \
    --calib_reference lambda_3.6kb.fasta \
    --hp_correct 1 \
    --trim_strategy 'dna' \
    --trim_threshold 2.5 \
    --trim_min_events 3 \
    --qscore_filtering \
    --min_qscore 7 \
    --num_callers 4 \
    --gpu_runners_per_device 2 \
    --chunk_size 1000 \
    --chunks_per_runner 1000 \
    --device "cuda:0" \
    --num_barcode_threads "$cpu" \
    --barcode_kits "$bk" \
    --trim_barcodes
    # --require_barcodes_both_ends

# If 1D2 experiment, run guppy_basecaller_1d2
guppy_basecaller_1d2     \
    --input_path "$basecalled"     \
    --save_path "$fastq"     \
    --index_file "$basecalled"/sequencing_summary.txt \
    --records_per_fastq 0 \
    --compress_fastq \
    --disable_pings \
    --config "${guppy_data_folder}"/dna_r9.5_450bps_1d2_raw.cfg \
    --num_callers 4 \
    --gpu_runners_per_device 2 \
    --chunk_size 1000 \
    --chunks_per_runner 1000  \
    --device "cuda:0"

# Trim and demultiplex already basecalled fastq
guppy_barcoder \
    --input_path "$basecalled" \
    --save_path "$fastq" \
    --recursive \
    --records_per_fastq 0 \
    --compress_fastq \
    --device "cuda:0" \
    --barcode_kits "$bk" \
    --trim_barcodes

# # no barcodes
# for i in $(find "$basecalled" -mindepth 1 -maxdepth 1 -type d); do
#     flag=$(basename "$i")
#     cat "${i}"/*.fastq.gz > "${i}"/"${flag}".fastq.gz
#     # Remove non compressed files
#     find "$i" -type f -name "*fastq_runid_*" -exec rm {} \;
# done

# Merge the basecalled and demultiplexed fastq
for i in $(find "$basecalled" -mindepth 2 -maxdepth 2 -type d); do  # pass and fail
    barcode=$(basename "$i")
    flag=$(basename $(dirname "$i"))
    # echo "${i}"/"${barcode}"_"${flag}".fastq.gz
    # cat "${i}"/*.fastq.gz > "${i}"/"${barcode}"_"${flag}".fastq.gz  #bash: /bin/cat: Argument list too long
    find "$i" -mindepth 1 -maxdepth 1 -type f -name "*.fastq.gz" -name "fastq_runid_*" \
        -exec cat {} \; > "${i}"/"${barcode}"_"${flag}".fastq.gz
    # # Remove non compressed files
    find "$i" -type f -name "*fastq_runid_*" -exec rm {} \;
done

#Rename samples
# Needs a translation files with 2 tab-separated columns
# barcode01 Sample1

# Parse translation talbe
#Do hash with conversion table
declare -A myArray=()

while IFS= read -r line || [[ -n "$line" ]]; do
    line="$(sed -e 's/ /\t/' -e 's/\r//g' -e 's/\n//g' <<< "$line")" #transform the space output field separator from read into tabs, remove carriage return
    key="$(cut -f 1 <<< "$line")"
    value="$(cut -f 2 <<< "$line")"
    myArray["${key}"]="${value}"
done < "$bc_desc"

# rename files
find "${basecalled}"/ -type f -name "*_pass*" -o -name "*_fail*" | while read i; do
    # echo "$i"
    pathPart="$(dirname "$i")"
    # echo "$pathPart"
    oldName="$(basename "$i")"
    # echo "$oldName"

    #for each file, check if a part of the name matches on
    for j in "${!myArray[@]}"; do
        # echo "$j"
        if [ "$(echo "$oldName" | grep "$j")" ]; then
            newName="$(echo "$oldName" | sed "s/"$j"/"${myArray["$j"]}"/")"
            fullNewName=""${pathPart}"/"${newName}""

            if [ -e "$rename" ]; then
                echo "Cannot rename "$oldName" to "$newName", file already exists. Skipping"
                continue
                # exit 1
            fi

            echo ""$i" -> "$fullNewName""
            mv "$i" "$fullNewName"
        fi
        if [ "$(echo "$pathPart" | grep "$j")" ]; then
            # rename folder too
            echo ""$pathPart" -> $(dirname "$pathPart")/"${myArray["$j"]}""
            mv $pathPart $(dirname "$pathPart")/"${myArray["$j"]}"
        fi
    done
done


##########
#        #
#   QC   #
#        #
##########


# # nanoQC on "sequencing_summary.txt" file
# [ -d "${qc}"/nanoQC/raw/summary ] || mkdir -p "${qc}"/nanoQC/raw/summary
# python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v3.1.5.py \
#     -s "${basecalled}"/sequencing_summary.txt \
#     -o "${qc}"/nanoQC/raw/summary

# Delete calibration strand folder if present
[ -d "${basecalled}"/calibration_strands ] && rm -rf "${basecalled}"/calibration_strands

# nanoQC on fastq files
[ -d "${qc}"/nanoQC/raw/fastq ] || mkdir -p "${qc}"/nanoQC/raw/fastq
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v3.1.5.py \
    -f "$basecalled" \
    -o "${qc}"/nanoQC/raw/fastq

function run_fastqc()
{
    # sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample=$(basename "${1%_filtered.fastq.gz}")
    sample="${sample%_pass.fastq.gz}"
    sample="${sample%_fail.fastq.gz}"

    [ -d "${2}"/"$sample" ] || mkdir -p "${2}"/"$sample"

    fastqc \
        --o "${2}"/"$sample" \
        --noextract \
        --threads $((cpu/maxProc)) \
        "$1"
}

export -f run_fastqc

# Create folder to store report
[ -d "${qc}"/fastqc/raw/pass ] || mkdir -p "${qc}"/fastqc/raw/pass
[ -d "${qc}"/fastqc/raw/fail ] || mkdir -p "${qc}"/fastqc/raw/fail

#raw
find "${basecalled}"/pass -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/raw/pass
find "${basecalled}"/fail -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/raw/fail

#Merge all FastQC reports together
source activate multiqc
multiqc \
    -o "${qc}"/fastqc/raw/pass \
    -n merged_pass.html \
    "${qc}"/fastqc/raw/pass
multiqc \
    -o "${qc}"/fastqc/raw/fail \
    -n merged_fail.html \
    "${qc}"/fastqc/raw/fail
source deactivate



################
#              #
#   Trimming   #
#              #
################


# Trim adapters
function chop()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    porechop \
        -i "$1" \
        -o "${trimmed}"/"${sample}"_trimmed.fastq.gz \
        --threads $((cpu/maxProc)) \
        | tee -a "${logs}"/trimming/"${sample}"_porechop.log

    # cleanup the log file by removing the progress lines
    cat "${logs}"/trimming/"${sample}"_porechop.log \
        | grep -vF '%' \
        > "${logs}"/trimming/"${sample}"_porechop.log.tmp

    mv "${logs}"/trimming/"${sample}"_porechop.log.tmp \
        "${logs}"/trimming/"${sample}"_porechop.log
}

export -f chop

[ -d "${logs}"/trimming ] || mkdir -p "${logs}"/trimming

find "${basecalled}/pass" -type f -name "*.fastq.gz" ! -name "*unclassified*" \
    | parallel  --bar \
                --env chop \
                --env trimmed \
                --env logs \
                --env cpu \
                --env maxProc \
                --jobs "$maxProc" \
                'chop {}'


#################
#               #
#   Filtering   #
#               #
#################


### Careful here. Filtered reads sometimes lead to a poorer quality assembly.
function filter()
{
    #sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample=$(basename "$1" "_trimmed.fastq.gz")

    # https://github.com/rrwick/Filtlong
    # --min_mean_q 90 -> phred 10
    # --min_length 1000 \
    # --target_bases -> max 100X coverage
    filtlong \
        --target_bases $((size*100)) \
        --keep_percent 90 \
        "$1" \
        | pigz > "${filtered}"/"${sample}"_filtered.fastq.gz
        # 2> >(tee "${logs}"/filtering/"${sample}".txt) \
}

export -f filter

[ -d "${logs}"/filtering ] || mkdir -p "${logs}"/filtering

find "$trimmed" -type f -name "*_trimmed.fastq.gz" \
    | parallel  --bar \
                --env filter \
                --env trimmed \
                --env size \
                --env filtered \
                --env logs \
                --jobs "$maxProc" \
                'filter {}'

# No filtering
# Just create symbolic links
# find "$basecalled" -type f -name "*_pass.fastq.gz" ! -name "*unclassified*" -exec ln -s {} "$filtered" \;


######################
#                    #
#   Quality Control  #
#                    #
######################


# nanoQC on fastq file
[ -d "${qc}"/nanoQC/filtered ] || mkdir -p "${qc}"/nanoQC/filtered
python3 /home/bioinfo/PycharmProjects/nanoQC/nanoQC_Guppy_v3.1.5.py \
    -f "$filtered" \
    -o "${qc}"/nanoQC/filtered

### FastQC###

# Create folder to store report
[ -d "${qc}"/fastqc/filtered ] || mkdir -p "${qc}"/fastqc/filtered


#filtered
find -L "$filtered" -type f -name "*.fastq.gz" \
    | parallel --bar run_fastqc {} "${qc}"/fastqc/filtered

#Merge all FastQC reports together
source activate multiqc
multiqc \
    -o "${qc}"/fastqc/filtered \
    -n merged_filtered.html \
    "${qc}"/fastqc/filtered
source deactivate


#########################################
#                                       #
#   Taxonomic Classification of reads   #
#                                       #
#########################################



[ -d "${qc}"/kraken2/raw ] || mkdir -p "${qc}"/kraken2/raw


function run_kraken2 ()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

    # Create folder to store report
    [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

    #run Kraken
    kraken2 \
        --memory-mapping \
        --db "$kraken2_db" \
        --output "${2}"/"${sample}"/"${sample}".kraken \
        --report "${2}"/"${sample}"/"${sample}".report.tsv \
        --threads "$cpu" \
        --gzip-compressed \
        "$1" \
        &> >(tee "${2}"/"${sample}"/"${sample}".kraken2.log)

    #Prepare result for display with Krona
    ktImportTaxonomy \
        -q 2 \
        -t 3 \
        -o "${2}"/"${sample}"/"${sample}".html \
        "${2}"/"${sample}"/"${sample}".kraken 
}

for i in $(find -L "${basecalled}/pass" -type f -name "*fastq.gz"); do
    run_kraken2 "$i" "${qc}"/kraken2/raw
done

# function run_centrifuge()
# {
#     sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

#     [ -d "${2}"/"${sample}" ] || mkdir -p "${2}"/"${sample}"

#     #build the command
    
#     centrifuge \
#         -p "$cpu" \
#         --mm \
#         -t \
#         --seed "$RANDOM" \
#         -x "$db" \
#         -U "$1" \
#         --report-file "${2}"/"${sample}"/"${sample}"_report.tsv \
#         > "${2}"/"${sample}"/"${sample}".tsv

#     cat "${2}"/"${sample}"/"${sample}".tsv | \
#         cut -f 1,3 | \
#         ktImportTaxonomy /dev/stdin -o "${2}"/"${sample}"/"${sample}".html

#     # Visualize the resutls in Firefow browser
#     firefox file://"${2}"/"${sample}"/"${sample}".html &
# }

# for i in $(find "${basecalled}"/pass -type f -name "*fastq.gz"); do
#     run_centrifuge "$i" "${qc}"/centrifuge/raw
# done


# ### CCMetagen
# # https://github.com/vrmarcelino/CCMetagen
# conda activate ccmetagen

# function ccmeta_it()
# {
#     input="$1"
#     sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
#     [ -d "${qc}"/ccmetagen/"$sample" ] || mkdir -p "${qc}"/ccmetagen/"$sample"


#     # Map reads with KMA
#     kma \
#         -i "$input" \
#         -o "${qc}"/ccmetagen/"${sample}"/"$sample" \
#         -t_db /media/30tb_raid10/db/CCMetagen/compress_ncbi_nt/ncbi_nt \
#         -t "$cpu" \
#         -bcNano \
#         -1t1 -mem_mode -and

#     CCMetagen.py \
#         -r nt \
#         -i "${qc}"/ccmetagen/"${sample}"/"${sample}".res \
#         -o "${qc}"/ccmetagen/"${sample}"/"${sample}".ccmeta

#     CCMetagen_merge.py \
#         -t Family \
#         -i "${qc}"/ccmetagen/"${sample}" \
#         -o "${qc}"/ccmetagen/"${sample}"/merged_samples
# }



# source deactivate


########################
#                      #
#   de novo assembly   #
#                      #
########################


function assemble()
{
    sample=$(basename "${1%_filtered.fastq.gz}")
    ass="$2"

    [ -d "${assemblies}"/"${ass}" ] || mkdir -p "${assemblies}"/"${ass}"
    [ -d "${qc}"/assembly_graphs/"$ass" ] || mkdir -p "${qc}"/assembly_graphs/"$ass"
    
    if [[ "$ass" == "unicycler" ]]; then
        
        unicycler \
            -l "$1" \
            -o "${assemblies}"/"${ass}"/"$sample" \
            -t $((cpu/maxProc)) \
            --verbosity 2 \
            --mode normal

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.gfa \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

        Bandage image \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
            "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

    elif [[ "$ass" == "unicycler_hybrid" ]]; then
        # $3 unmerged R1
        # $4 unmerged R2
        # $5 merged
        unicycler \
            -1 "$3" \
            -2 "$4" \
            -s "$5" \
            --no_correct \
            -l "$1" \
            -o "${assemblies}"/"${ass}"/"$sample" \
            -t $((cpu/maxProc)) \
            --verbosity 2 \
            --mode normal \
            --pilon_path "${prog}"/pilon/pilon-dev.jar

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta
    elif [[ "$ass" == "canu" ]]; then
        canu \
            -p "$sample" \
            -d "${assemblies}"/"${ass}"/"$sample" \
            genomeSize="$size" \
            maxThreads=$((cpu/maxProc)) \
            -nanopore-raw "$1"

        # rename assembly
        mv "${assemblies}"/"${ass}"/"${sample}"/"${sample}".contigs.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

        mv "${assemblies}"/"${ass}"/"${sample}"/"${sample}".contigs.gfa \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

        Bandage image \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
            "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

        # Cleanup
        find "${assemblies}"/"${ass}"/"$sample" -mindepth 1 -maxdepth 1 \
            ! -name "*.fasta" \
            ! -name "*.gfa" \
            ! -name "*.report" \
            -o -name "*unitigs*" \
            -o -name "*unassembled*" \
            -exec rm -rf {} \;

    elif [[ "$ass" == "flye" ]]; then
        flye \
            --nano-raw "$1" \
            --out-dir "${assemblies}"/"${ass}"/"$sample" \
            --genome-size "$size" \
            --iterations 3 \
            --threads $((cpu/maxProc))

        # Rename assembly
        mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

        mv "${assemblies}"/"${ass}"/"${sample}"/assembly_graph.gfa \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

        Bandage image \
            "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
            "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

        # Maybe some cleanup
        find "${assemblies}"/"${ass}"/"$sample" -mindepth 1 -maxdepth 1 \
            ! -name "*.fasta" \
            ! -name "*.gfa" \
            ! -name "*.log" \
            ! -name "*.txt" \
            -exec rm -rf {} \;
    # else
    #     echo 'Please chose one of the following assembler: "unicycler", "unicycler_hybrid", "canu" or "flye"'
    fi
}

export -f assemble

find -L "${filtered}" -type f -name "*.fastq.gz" \
    | parallel  --bar \
                --env cpu \
                --env maxProc \
                --env assemble \
                --env assembler \
                --env assemblies \
                --env fastq2fasta \
                --env size \
                --jobs "$maxProc" \
                "assemble {} "$assembler""


function fastq2fasta ()
{
    input_fastq="$1"  # assume gzipped file

    if [ $(echo "$input_fastq" | grep -F ".fastq.gz") ]; then
        zcat "$input_fastq" \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > "${input_fastq%.fastq.gz}.fasta"
    elif [ $(echo "$input_fastq" | grep -F ".fastq") ]; then
        cat "$input_fastq" \
        | sed -n '1~4s/^@/>/p;2~4p' \
        > "${input_fastq%.fastq}.fasta"
    else
        echo "Wrong file extension. Please use \".fastq\" or \".fastq.gz\""
        exit 1
    fi
}

export -f fastq2fasta


function shasta_it()
{
    sample=$(basename "${1%_filtered.fastq.gz}")
    ass=shasta

    file_ext=""

    if [ $(echo "$1" | grep -F ".fastq.gz") ]; then
        file_ext=".fastq.gz"
    elif [ $(echo "$1" | grep -F ".fastq") ]; then
        file_ext=".fastq"
    else
        echo "Wrong file extension. Please use \".fastq\" or \".fastq.gz\""
    # exit 1
    fi

    # Convert reads to fasta
    echo "Converting fastq to fasta for Shasta..."
    fastq2fasta "$1"
    [ -d "${assemblies}"/"${ass}" ] || mkdir -p "${assemblies}"/"${ass}"
    mv "${1%$file_ext}.fasta" "${assemblies}"/"${ass}"/"${sample}".fasta

    # Two problems here:
    #   1- Requires sudo
    #   2- Defaults to all CPUs
    # shasta --command assemble \
    #     --input "${assemblies}"/"${ass}"/"${sample}".fasta \
    #     --assemblyDirectory "${assemblies}"/"${ass}"/"$sample" \
    #     --memoryBacking 2M \
    #     --memoryMode filesystem \
    #     --Reads.minReadLength=10000
    shasta \
        --command assemble \
        --input ""${assemblies}"/"${ass}"/"${sample}".fasta" \
        --assemblyDirectory "${assemblies}"/"${ass}"/"$sample" \
        --memoryBacking 2M \
        --memoryMode filesystem \
        --Reads.minReadLength=10000 \
        --threads "$cpu" \
        --MinHash.minBucketSize 5 \
        --MinHash.maxBucketSize 30 \
        --MinHash.minFrequency 5 \
        --Align.minAlignedFraction 0.4 \
        --Assembly.consensusCaller 'Bayesian:guppy-3.0.5-a'

    shasta --command cleanupBinaryData \
        --assemblyDirectory "${assemblies}"/"${ass}"/"$sample"

    # remove input fasta file
    rm "${assemblies}"/"${ass}"/"${sample}".fasta

    # Rename outputs
    mv "${assemblies}"/"${ass}"/"${sample}"/Assembly.fasta \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

    mv "${assemblies}"/"${ass}"/"${sample}"/Assembly.gfa \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

    # Make assembly graph
    [ -d "${qc}"/assembly_graphs/"$ass" ] || mkdir -p "${qc}"/assembly_graphs/"$ass"
    Bandage image \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa \
        "${qc}"/assembly_graphs/"${ass}"/"${sample}".png

    # Reformat fasta
    perl "${scripts}"/formatFasta.pl \
        -i "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta \
        -o "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta1
    mv "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta1 \
        "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta
}

# Loop all the samples
if [[ "$assembler" == "shasta" ]]; then
    for i in $(find -L "${filtered}" -type f -name "*.fastq.gz"); do
        shasta_it "$i" # Needs sudo
    done
fi


#################
#               #
#   Polishing   #
#               #
#################

# Marginpolish + HELEN to use on shasta assemblies -> seems to be more for human genome assemblies
# https://github.com/UCSC-nanopore-cgl/marginPolish
# https://github.com/kishwarshafin/helen


# Do 4 rounds of Racon polishing


# Polish with Medaka
# https://github.com/nanoporetech/medaka
source activate medaka
# source ~/my_virtualenv/medaka_gpu/bin/activate

function polish_medaka()
{
    sample=$(basename "$1" ".fasta")

    [ -d "${polished}"/"${assembler}"/"$sample" ] || mkdir -p "${polished}"/"${assembler}"/"$sample"

    # It is crucially important to specify the correct model (-m) according to the basecaller used.
    # Allowed values can be found by running "medaka tools list\_models"
    # -i "${basecalled}"/pass/"${sample}"/"${sample}"_pass.fastq.gz \
    medaka_consensus \
        -i "${filtered}"/"${sample}"_filtered.fastq.gz \
        -d "$1" \
        -o "${polished}"/"${assembler}"/"$sample" \
        -t $((cpu/maxProc)) \
        -m r941_min_high_g344 \
        2>&1 | tee "${polished}"/"${assembler}"/"${sample}"/medaka.log

    # Rename sample
    mv "${polished}"/"${assembler}"/"${sample}"/consensus.fasta \
        "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta

    # Remove temporary files
    find "${polished}"/"${assembler}"/"${sample}" -type f \
        ! -name "*_medaka.fasta" \
        ! -name "*.gfa" \
        ! -name "*.log" \
        -exec rm {} \;

    # Reformat fasta
    perl "${scripts}"/formatFasta.pl \
        -i "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta \
        -o "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta1
    mv "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta1 \
        "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta
}

export -f polish_medaka

find "${assemblies}"/"$assembler" -maxdepth 2 -type f -name "*.fasta" \
    | parallel  --bar \
                --env polish_medaka \
                --env polished \
                --env cpu \
                --env assembler \
                --env maxProc \
                'polish_medaka {}'

source deactivate medaka
# deactivate

# Compare pre- and post-nanopolished genomes
function compare_assemblies()
{
    # sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample=$(basename "$1" "_medaka.fasta")

    [ -d "${qc}"/mummer/"${assembler}"/"$sample" ] || mkdir -p "${qc}"/mummer/"${assembler}"/"$sample"

    nucmer \
        --mum \
        --threads=$((cpu/maxProc)) \
        --delta "${qc}"/mummer/"${assembler}"/"${sample}"/"${sample}".delta \
        "${polished}"/"${assembler}"/"${sample}"/"${sample}"_medaka.fasta \
        "${assemblies}"/"${assembler}"/"${sample}"/"${sample}".fasta

    mummerplot \
        -p "${qc}"/mummer/"${assembler}"/"${sample}"/"${sample}" \
        -title 'Pre versus Post medaka' \
        --layout \
        --large \
        --color \
        --png \
        "${qc}"/mummer/"${assembler}"/"${sample}"/"${sample}".delta

    dnadiff \
        -p "${qc}"/mummer/"${assembler}"/"${sample}"/"${sample}" \
        -d "${qc}"/mummer/"${assembler}"/"${sample}"/"${sample}".delta
}

export -f compare_assemblies

find "${polished}"/"$assembler" -type f -name "*.fasta" | \
parallel    --bar \
            --env compare_assemblies \
            --env polished \
            --env assemblies \
            --env assembler \
            --env qc \
            --env cpu \
            --env maxProc \
            --jobs "$maxProc" \
            "compare_assemblies {}"


# Fix assembly start -> put dnaA gene first
function fix_start()
{
    sample=$(basename "$1" "_medaka.fasta")
    [ -d "${baseDir}"/fixstart/"$assembler" ] || mkdir -p "${baseDir}"/fixstart/"$assembler"

    circlator fixstart \
        "$1" \
        "${baseDir}"/fixstart/"$assembler"/"$sample"

    # check if dnaA was found
    # if not, dont keep the new "fixed" assebmly because it has be rerooted from the middle


    # cleanup
    find "${baseDir}"/fixstart/"$assembler" -type f ! -name "*.fasta" -exec rm {} \;
}

export -f fix_start

find "${polished}"/"$assembler" -type f -name "*.fasta" \
| parallel --bar --env fix_start --env baseDir 'fix_start {}'




### Coverage ###
function get_coverage()  # unsing unmerged reads only
{
    # sample=$(basename "$1" | cut -d '_' -f 1)
    sample=$(basename "$1" "_medaka.fasta")
    fastq_file="${filtered}"/"${sample}"_filtered.fastq.gz

    [ -d "${qc}"/coverage/"${assembler}"/"$sample" ] || mkdir -p "${qc}"/coverage/"${assembler}"/"$sample"

    minimap2 \
        -ax map-ont \
        -t $((cpu/maxProc)) \
        "$1" \
        "$fastq_file" | \
    samtools view -@ $((cpu/maxProc)) -b -h -F 4 - | \
    samtools sort -@ $((cpu/maxProc)) - | \
    samtools rmdup - "${qc}"/coverage/"${assembler}"/"${sample}"/"${sample}".bam

    # Index bam file
    samtools index -@ $((cpu/maxProc)) "${qc}"/coverage/"${assembler}"/"${sample}"/"${sample}".bam

    #Average genome depth of coverage
    average_cov=$(samtools depth \
        "${qc}"/coverage/"${assembler}"/"${sample}"/"${sample}".bam  \
        | awk '{sum+=$3} END { print sum/NR}')

    printf "%s\t%.*f\n" "$sample" 0 "$average_cov" | tee -a "${qc}"/coverage/"${assembler}"/average_cov.tsv

    #Remove reads in fasta format  -> for nanopolish
    # rm "${polished}"/"${sample}"/"${sample}"_reads.fasta
}

export -f get_coverage

[ -d "${qc}"/coverage/"${assembler}" ] || mkdir -p "${qc}"/coverage/"${assembler}"
echo -e "Sample\tAverage_Cov" > "${qc}"/coverage/"${assembler}"/average_cov.tsv

find "${polished}"/"$assembler" -type f -name "*.fasta" | \
    parallel    --bar \
                --env get_coverage \
                --env cpu \
                --env maxProc \
                --env assembler \
                --env qc \
                --jobs "$maxProc" \
                "get_coverage {}"

# Qualimap
function run_qualimap()
{
    sample=$(basename "$1" ".bam")
    
    [ -d "${qc}"/coverage/"${assembler}"/qualimap/"$sample" ] || mkdir -p "${qc}"/coverage/"${assembler}"/qualimap/"$sample"

    qualimap bamqc \
        --paint-chromosome-limits \
        -bam "$1" \
        --java-mem-size="${mem}"G \
        -nt $((cpu/maxProc)) \
        -outdir "${qc}"/coverage/"${assembler}"/qualimap/"$sample" \
        -outfile "${sample}" \
        -outformat HTML

    # Remove bam files
    rm -rf "${qc}"/coverage/"${assembler}"/"$sample"
}

export -f run_qualimap

find "${qc}"/coverage/"${assembler}" -type f -name "*.bam" |
parallel    --bar \
            --env run_qualimap \
            --env qc \
            --env mem \
            --env cpu \
            --env assembler \
            --env maxProc \
            --jobs "$maxProc" \
            'run_qualimap {}'


# Blast genomes on nr
function blast()
{
    sample=$(basename "$1" "_medaka.fasta")

    [ -d "${qc}"/blast/"${assembler}"/"$sample" ] || mkdir -p "${qc}"/blast/"${assembler}"/"$sample"

    blastn \
        -db nt \
        -query "$1" \
        -out "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv \
        -evalue "1e-10" \
        -outfmt '6 qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames sskingdoms staxids' \
        -num_threads $((cpu/maxProc)) \
        -max_target_seqs 20

    # Add header
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv.tmp

    # Sort hits
    cat "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv \
        | sed '1d' \
        | sort -t $'\t' -k1,1g -k13,13gr \
        >> "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv.tmp

    # Replace original file
    mv "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv.tmp \
        "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv

    # Best hits only
    echo -e "qseqid\tsseqid\tstitle\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tsscinames\tsskingdoms\tstaxids" \
        > "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp

    # Only works with Unicycler
    # TODO -> Change so it's more universal (works with other assembler headers)
    cat "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.all.tsv \
        | sed '1d' \
        | sort -t $'\t' -k1,1g -k13,13gr \
        | sort -t $'\t' -uk1,1g \
        >> "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp

    mv "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.bestHit.tsv.tmp \
        "${qc}"/blast/"${assembler}"/"${sample}"/"${sample}".blastn.bestHit.tsv
}

# export -f blast

# find "$polished" -type f -name "*_nanopolished.fasta" \
#     | parallel  --bar \
#                 --env blast \
#                 --env cpu \
#                 --env maxProc \
#                 --env qc \
#                 --jobs "$maxProc" \
#                 "blast {}"

# It's faster to run the blast on long contigs with all cores one sample at the time
# than to parallel blast the samples with fewer cores, at least for nt
c=0
n=$(find "${polished}"/"$assembler" -type f -name "*.fasta" | wc -l)
for i in $(find "${polished}"/"$assembler" -type f -name "*.fasta"); do
    let c+=1
    echo -ne "Blasting assembly "${c}"/"${n}" \\r"
    blast "$i"
done


### QUAST ###

declare -a genomes=()
for i in $(find "${polished}"/"$assembler" -type f -name "*.fasta"); do 
    genomes+=("$i")
done

source activate quast

quast.py \
    --output-dir "${qc}"/quast/"${assembler}"/all \
    --threads "$cpu" \
    --min-contig "$smallest_contig" \
    --est-ref-size "$size" \
    ${genomes[@]}


# Make quast report on individual assembly
function run_quast()
{
    sample=$(basename "$1" "_medaka.fasta")

    quast.py \
        -m "$smallest_contig" \
        -t $((cpu/maxProc)) \
        -o "${qc}"/quast/"${assembler}"/"$sample" \
        --min-contig "$smallest_contig" \
        --est-ref-size "$size" \
        $1  # don't put in quotes
}

#make function available to parallel
export -f run_quast  # -f is to export functions

#run paired-end merging on multiple samples in parallel
find "${polished}"/"$assembler" -type f -name "*.fasta" \
    | parallel  --bar \
                --env run_quast \
                --env maxProc \
                --env cpu \
                --env qc \
                --env assembler \
                --env smallest_contig \
                --jobs "$maxProc" \
                "run_quast {}"

source deactivate


###################
#                 #
#   Annotation    #
#                 #
###################


source activate prokka

function annotate()
{
    sample=$(basename "$1" "_medaka.fasta")
    sample="${sample%.fasta}"

    #Prokka
    prokka  --outdir "${annotation}"/"${assembler}"/"$sample" \
            --force \
            --prefix "$sample" \
            --kingdom "$kingdom" \
            --genus "$genus" \
            --species "$species" \
            --strain "$sample" \
            --gram "$gram" \
            --locustag "$locustag" \
            --compliant \
            --centre "$centre" \
            --cpus $((cpu/maxProc)) \
            --rfam \
            --protein /media/30tb_raid10/db/blast/prokka_hypoth/Psyringae_tomato_hypth.fasta \
            "$1"

    #extract hypothetical proteins
    cat "${annotation}"/"${assembler}"/"${sample}"/"${sample}".faa | \
        awk '{if(substr($0,1,1)==">"){if (p){print "\n";} print $0} else printf("%s",$0);p++;}END{print "\n"}' | \
        grep --no-group-separator -A 1 -F "hypothetical protein" \
        > "${annotation}"/"${assembler}"/"${sample}"/"${sample}"_hypoth.faa

    echo -e ""$sample" hypothetical proteins (round1): $(cat "${annotation}"/"${assembler}"/"${sample}"/"${sample}".faa | grep -ic "hypothetical")" \
        | tee -a "${logs}"/log.txt
}

export -f annotate

find "${polished}"/"$assembler" -type f -name "*.fasta" | \
    parallel    --bar \
                --env annotate \
                --env annotation \
                --env kingdom \
                --env genus \
                --env species \
                --env gram \
                --env locustag \
                --env centre \
                --env cpu \
                --env maxProc \
                --env assembler \
                --env scripts \
                --jobs "$maxProc" \
                "annotate {}"

source deactivate

# Run all the assemblers in a row
# for i in unicycler shasta flye canu; do
#     export assembler="$i"
# done


### Resfinder

function run_resfinder ()
{
    # https://bitbucket.org/genomicepidemiology/resfinder/overview
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample="${sample%.fasta}"

    [ -d "${amr}"/resfinder/"${assembler}"/"$sample" ] || mkdir -p "${amr}"/resfinder/"${assembler}"/"$sample"

    python3 "${prog}"/resfinder/resfinder.py \
        -i "$1" \
        -o "${amr}"/resfinder/"${assembler}"/"$sample" \
        -p "${prog}"/resfinder/resfinder_db/ \
        -t 0.9 \
        -l 0.6 #\
        # 1> >(tee "${amr}"/resfinder/"${sample}"/"${sample}"_resfinder.txt)

    rm -rf "${amr}"/resfinder/"${assembler}"/"${sample}"/tmp
}

export -f run_resfinder

find "${polished}"/"$assembler" -type f -name "*.fasta" | \
    parallel    --env run_resfinder \
                --env resfinder_db \
                --env assembler \
                "run_resfinder {}"

# Create merged report
echo -e 'Sample\tResistance gene\tIdentity\tAlignment Length/Gene Length\tCoverage\tPosition in reference\tContig\tPosition in contig\tPhenotype\tAccession no.' \
    > "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv.tmp
    
for i in $(find "${amr}"/resfinder/"$assembler" -name "*results_tab.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))
    sample="${sample%.fasta}"

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv.tmp
done

# sort by sample name (column 1), then by Identity (column 3)
(cat "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv.tmp | head -n 1;
    cat "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k3,3) \
    > "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv

rm "${amr}"/resfinder/"${assembler}"/resfinder_merged.tsv.tmp



### RGI (CARD)

function run_rgi()
{
    sample=$(cut -d "_" -f 1 <<< $(basename "$1"))
    sample="${sample%.fasta}"

    [ -d "${amr}"/rgi/"${assembler}"/"$sample" ] || mkdir -p "${amr}"/rgi/"${assembler}"/"$sample"

    "${prog}"/rgi-4.2.2/./rgi main \
        -i "$1" \
        -o  "${amr}"/rgi/"${assembler}"/"${sample}"/"$sample" \
        -t 'contig' \
        -a 'BLAST' \
        -n $((cpu/maxProc)) \
        --clean \
        -d "chromosome"
}

export -f run_rgi

find "${polished}"/"$assembler" -type f -name "*.fasta" \
    | parallel    --bar \
                --env run_rgi \
                --env resfinder_db \
                --env assembler \
                "run_rgi {}"

# Create merged report
echo -e 'Sample\tORF_ID\tContig\tStart\tStop\tOrientation\tCut_Off\tPass_Bitscore\tBest_Hit_Bitscore\tBest_Hit_ARO\tBest_Identities\tARO\tModel_typeSNPs_in_Best_Hit_ARO\tOther_SNPs\tDrug\tClass\tResistance\tMechanism\tAMR\tGene\tFamily\tPredicted_DNA\tPredicted_Protein\tCARD_Protein_Sequence\tPercentage_Length_of_Reference_Sequence\tID\tModel_ID' \
    > "${amr}"/rgi/"${assembler}"/rgi_merged.tsv.tmp
    
for i in $(find "${amr}"/rgi/"${assembler}" -name "*.txt"); do
    # sample name is folder name
    sample=$(basename $(dirname "$i"))
    sample="${sample%.fasta}"

    # Add a leading column with sample name
    cat "$i" \
        | sed -e '1d' \
        | awk -F $'\t' -v s="$sample" 'BEGIN {OFS = FS} {print s,$0}' \
        >> "${amr}"/rgi/"${assembler}"/rgi_merged.tsv.tmp
done

# sort by Sample name (column 1), then by Cutt_off (column 7)
(cat "${amr}"/rgi/"${assembler}"/rgi_merged.tsv.tmp | head -n 1;
    cat "${amr}"/rgi/"${assembler}"/rgi_merged.tsv.tmp | sed -e '1d' | sort -t $'\t' -k1,1 -k7,7) \
    > "${amr}"/rgi/"${assembler}"/rgi_merged.tsv

rm "${amr}"/rgi/"${assembler}"/rgi_merged.tsv.tmp


### Phaster

#trim assemblies
function phaster_trim()
{
    sample=$(cut -d '_' -f 1 <<< $(basename "$1"))
    sample="${sample%.fasta}"

    # http://phaster.ca/instructions
    if [ $(cat "$1" | grep -Ec "^>") -gt 1 ]; then  # if more than one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            2000 \
            "$1" \
            > "${phaster}"/"${assembler}"/assemblies/"${sample}".fasta
    elif [ $(cat "$1" | grep -Ec "^>") -eq 1 ]; then  # if only one contig
        #remove contigs smaller than 2000 bp from assembly
        perl "${prog}"/phage_typing/removesmallscontigs.pl \
            1500 \
            "$1" \
            > "${phaster}"/"${assembler}"/assemblies/"${sample}".fasta
    else
        echo "No assembly for "$sample""  # Should not get here!
        # exit 1
    fi
}

#make function available to parallel
export -f phaster_trim  # -f is to export functions

 # To store trimmed assemblies for phaster submission
[ -d "${phaster}"/"${assembler}"/assemblies ] || mkdir -p "${phaster}"/"${assembler}"/assemblies 

#run trimming on multiple assemblies in parallel
find "${polished}"/"$assembler" -type f -name "*.fasta" \
    | parallel  --bar \
                --env phaster_trim \
                --env prog \
                'phaster_trim {}'

# Get phaster results
# python3 ~/scripts/checkPhasterServer.py --submit --check \
#     -i "${phaster}"/assemblies \
#     -o "$phaster"

# Run local version of phaster
# /media/30tb_raid10/db/blast/prophages/prophage_virus_filtered_fixedheaders.db
for i in $(find "${phaster}"/"${assembler}"/assemblies -type f -name "*.fasta"); do
    sample=$(basename "$1" ".fasta")
    [ -d "${phaster}"/"${assembler}"/"$sample" ] || mkdir -p "${phaster}"/"${assembler}"/"$sample"

    perl /home/bioinfo/prog/phaster-app/scripts/phaster.pl \
        -c \
        -i "$i" \
        -o "${phaster}"/"${assembler}"/"${sample}"/"$sample"
done

# Cleanup
rm -rf "${phaster}"/"${assembler}"/assemblies


#reformat phaster output
function reformat_output ()
{
    # #unzip phaster result file
    # [ -d "${phaster}"/"$sample" ] || mkdir -p "${phaster}"/"$sample"
    # unzip "$1" -d "${phaster}"/"${sample}"

    #Convert summary.txt to tsv format
    cat "${1}"/summary.txt \
        | sed '1,/gi|/d' \
        | awk '{$1=$1;print}' \
        | sed '/^-/d' \
        | tr " " "\t" \
        > "${1}"/summary.tsv

    #Convert detail.txt to tsv format
    cat "${1}"/detail.txt \
        | sed -e '/^-/d' -e 's/  /@/g' \
        | tr -s '@' \
        | sed 's/@ /@/' \
        | tr -s '[:space:]+' \
        | tr "@" "\t" \
        > "${1}"/detail.tsv
}

for i in $(find "${phaster}"/"$assembler" -mindepth 1 -type d ! -name assemblies); do
    reformat_output "$i"
done

# Rename output file
for i in $(find "${phaster}"/"$assembler" -type f -name "region_DNA.txt"); do
    d=$(dirname "$i")
    sample=$(basename "$d")

    mv "$i" "${d}"/"${sample}"_phages.fasta
done

#Add sample name and entry number to fasta header
for i in $(find "${phaster}"/"$assembler" -name "*_phages.fasta"); do
    sed -i "s/^>/>"${sample}"_/" "$i"
done












# #######
# # assemble metagenomic data

# function meta_assemble()
# {
#     ass=flye_meta
#     sample=$(cut -d "_" -f 1 <<< $(basename "$1"))

#     [ -d "${assemblies}"/flye_meta/"$sample" ] || mkdir -p "${assemblies}"/flye_meta/"$sample"

#     flye \
#         --meta \
#         --nano-raw "$1" \
#         --out-dir "${assemblies}"/"${ass}"/"$sample" \
#         --genome-size "$size" \
#         --iterations 3 \
#         --threads $((cpu/maxProc))

#     # Rename assembly
#     mv "${assemblies}"/"${ass}"/"${sample}"/assembly.fasta \
#         "${assemblies}"/"${ass}"/"${sample}"/"${sample}".fasta

#     mv "${assemblies}"/"${ass}"/"${sample}"/assembly_graph.gfa \
#         "${assemblies}"/"${ass}"/"${sample}"/"${sample}".gfa

#     # Maybe some cleanup
#     find "${assemblies}"/"${ass}"/"$sample" -mindepth 1 -maxdepth 1 \
#         ! -name "*.fasta" \
#         ! -name "*.gfa" \
#         ! -name "*.log" \
#         ! -name "*.report" \
#         -o -name "*unitigs*" \
#         -exec rm -rf {} \;
# }

# export -f meta_assemble

# find "${basecalled}"/pass -type f -name "*.fastq.gz" \
#     | parallel  --bar \
#                 --env cpu \
#                 --env maxProc \
#                 --env meta_assemble \
#                 --env assembler \
#                 --env assemblies \
#                 --env size \
#                 --jobs "$maxProc" \
#                 'meta_assemble {}'


# # Mapping to reference
# # end to end
# function map_bowtie2()
# {
#     sample=$(cut -d '_' -f 1 <<< $(basename "$2"))

#     bowtie2 \
#         --mm \
#         --very-sensitive \
#         --threads $((cpu/maxProc)) \
#         -x "$1" \
#         -U "$2" \
#     | samtools view -bh -F4 -@ $((cpu/maxProc)) - \
#     | samtools sort -@ $((cpu/maxProc)) - \
#     > "${baseDir}"/mapped/"${sample}".bam

#     samtools index "${baseDir}"/mapped/"${sample}".bam
# }

# export -f map_bowtie2

# find "${basecalled}"/pass -name "*.fastq.gz" \
# | parallel  --bar \
#             --env map_bowtie2 \
#             --env cpu \
#             --env maxProc \
#             map_bowtie2 /media/30tb_raid10/ref/crypto/Cryptosporidium_parvum_strain_Iowa_II {}


