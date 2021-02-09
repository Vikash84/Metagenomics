files=(
Bird18
Bird19
Bird20
Bird21
Bird22

)

for file in "${files[@]}";
    do
	./Trinity --seqType fq --max_memory 220G --left /home/mang/project/bird/"$file"_L001_R1.fastq.gz,/home/mang/project/bird/"$file"_L002_R1.fastq.gz --right /home/mang/project/bird/"$file"_L001_R2.fastq.gz,/home/mang/project/bird/"$file"_L002_R2.fastq.gz --CPU 12 --full_cleanup --output "$file".trinity --normalize_reads --trimmomatic    ###Trinity

	./util/align_and_estimate_abundance.pl --transcripts "$file".trinity.Trinity.fasta --seqType fq --single /home/mang/project/bird/"$file"_L001_R1.fastq.gz,/home/mang/project/bird/"$file"_L002_R1.fastq.gz,/home/mang/project/bird/"$file"_L001_R2.fastq.gz,/home/mang/project/bird/"$file"_L002_R2.fastq.gz --est_method RSEM --aln_method bowtie2 --trinity_mode --output_dir "$file"_RSEM --thread_count 12 --prep_reference   ### RSEM abundance estimation
	
	sed "s/ path.*//gi" "$file".trinity.Trinity.fasta | sed 's/|/_/' | sed "s/ /_/gi" | sed "s/=//gi" | sed "s/TRINITY/"$file"/gi" > temp.fasta   ### Modify the name
	awk '{if (substr($0,1,1)==">"){print "\n"$0} else printf("%s",$0);p++;}END{print "\n"}' temp.fasta | sed 1d > "$file".trinity.fasta           ### Make oneline sequence
 
	blastn -query "$file".trinity.fasta -db /databases/ncbi/public/nt/nt -out "$file".trinity.fasta.nt -evalue 1E-10 -outfmt "6 qseqid qlen sacc salltitles pident length evalue sskingdoms" -max_target_seqs 3 -num_threads 16   ### blastn against nt database    #### Blastn
	diamond blastx -q "$file".trinity.fasta -d /databases/diamond/nr2 -o "$file".trinity.fasta.nr -e 1E-5 -k 3 -p 12 -f 6 qseqid qlen sseqid stitle pident length evalue --more-sensitive  #### Diamond blastx
    done;
