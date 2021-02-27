#!/bin/sh

#  ReadStats.sh
#  
#
#  Created by ellisrichardj on 10/10/2019.
#  

#  Define inputs - sample name

pair_id=$1

# Count reads in each catagory, using '+' as the only character on a line (^+$) as a proxy in fastq files 
# In fastq format '+' appears on the line between basecalls and quality scores for a given read

    raw_R1=$(zgrep -c "^+$" ${pair_id}_*_R1_*.fastq.gz) # counts number of reads in file
    rm ${pair_id}_*_R1_*.fastq.gz
    rm ${pair_id}_*_R2_*.fastq.gz
    uniq_R1=$(grep -c "^+$" ${pair_id}_uniq_R1.fastq) # counts number of reads in file
    rm `readlink ${pair_id}_uniq_R2.fastq`
    trim_R1=$(grep -c "^+$" ${pair_id}_trim_R1.fastq) # counts number of reads in file
    num_map=$(samtools view -c ${pair_id}.mapped.sorted.bam) # samtools counts the number of mapped reads
    samtools depth -a ${pair_id}.mapped.sorted.bam > depth.txt # samtools outputs depth at each position
    avg_depth=$(awk '{sum+=$3} END { print sum/NR}' depth.txt)
    zero_cov=$(awk '$3<1 {++count} END {print count}' depth.txt)
    sites=$(awk '{++count} END {print count}' depth.txt)
    rm depth.txt
    
# Caluclate values and percentages

    num_raw=$(($raw_R1*2))
    num_uniq=$(($uniq_R1*2))
    num_trim=$(($trim_R1*2))
    pc_aft_dedup=$(echo "scale=2; ($num_uniq*100/$num_raw)" |bc)
    pc_aft_trim=$(echo "scale=2; ($num_trim*100/$num_raw)" |bc)
    pc_mapped=$(echo "scale=2; ($num_map*100/$num_trim)" |bc)
    genome_cov=$(echo "scale=2; (100-($zero_cov*100/$sites))" |bc)

# Define thresholds for flag assignment

    mindepth=10 # require an average of at least 10 reads per site 
    minpc=60 # require at least 60% of data maps to genome
    minreads=600000 # require at least 600,000 raw reads per sample
    
# This section assigns 'flags' based on the number of reads and the proportion mapping to reference genome
    
    if [ ${avg_depth%%.*} -ge $mindepth ] && [ ${pc_mapped%%.*} -gt $minpc ]; then flag="Pass"
        elif [ ${avg_depth%%.*} -lt $mindepth ] && [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="Comtaminated"
        elif [ ${avg_depth%%.*} -lt $mindepth ] && [ $num_trim -lt $minreads ]; then flag="InsufficientData"
#        elif [ ${pc_mapped%%.*} -lt $minpc ] && [ $num_trim -gt $minreads ]; then flag="q_OtherMycobact"
        else flag="CheckRequired"
    fi

# Write values to csv file

    echo "Sample,NumRawReads,NumDedupReads,%afterDedup,NumTrimReads,%afterTrim,NumMappedReads,%Mapped,MeanDepth,GenomeCov,Outcome" > ${pair_id}_stats.csv
    echo "${pair_id},"$num_raw","$num_uniq","$pc_aft_dedup","$num_trim","$pc_aft_trim","$num_map","$pc_mapped","$avg_depth","$genome_cov","$flag"" >> ${pair_id}_stats.csv
    echo "$flag" > outcome.txt
