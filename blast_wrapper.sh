#!/bin/bash

Usage="""
Usage: blast_to_gff_wrapper.sh -q <query file> -d <database file> -p <blast program> -o <output file>
Options:
        -h      :       Help. What you are reading now.
        -q      :       Query. Put the path to your query fasta here.
        -d      :       Database. Put the path to your blast database here.
        -o      :       Output. Put the name or path and name to your output location here. Default: blastn
        -p      :       Program. Currently only confirmed to work with tblastn
                        but others should work too.
        -t      :       Threads. Number of threads/processors you want the blast analysis to
                        run on. (Default: 1)
"""
OUTPUT=blast.out
THREADS=2
PROGRAM=tblastn

while getopts :q:d:o:p:t:h opt; do
  case $opt in
        q)
                echo "-q (query) was input as $OPTARG" >&2
                QUERY=$OPTARG
        ;;
        d)
                echo "-d (database) was input as $OPTARG" >&2
                DATABASE=$OPTARG
        ;;
        o)
                echo "-o (output) was input as $OPTARG" >&2
                OUTPUT=$OPTARG
        ;;
        p)
                echo "-p (program) was input as $OPTARG" >&2
                PROGRAM=$OPTARG
        ;;
        t)
                echo "-t (threads) was input as $OPTARG" >&2
                THREADS=$OPTARG
        ;;
        l)
                echo "-l (long) was triggered, long output triggered" >&2
                LONG=true
        ;;
        k)
        echo "-k (keep) was triggered, blast output will be kept" >&2
        KEEP=true
        ;;
        h)
                echo "$Usage"
                exit 1
        ;;
        \?)
                echo "Invalid option: -$OPTARG" >&2
                echo "Type $0 -h for usage"
                exit 1
        ;;
  esac
done

if [[ $PROGRAM = "tblastn" ]] || [[ $PROGRAM = "TBLASTN" ]] && [[ $LONG = false ]] ; then
        makeblastdb -dbtype prot -in $QUERY -out $QUERY
        tblastn -query $QUERY -db $DATABASE -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" -out $OUTPUT -num_threads $THREADS


elif [[ $PROGRAM = "blastx" ]] || [[ $PROGRAM = "BLASTX" ]] && [[ $LONG = false ]] ; then
        makeblastdb -dbtype prot -in $QUERY -out $QUERY
        blastx -query $QUERY -db $DATABASE -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" -out $OUTPUT -num_threads $THREADS

elif [[ $PROGRAM = "blastn" ]] || [[ $PROGRAM = "BLASTN" ]] && [[ $LONG = false ]] ; then
        makeblastdb -dbtype nucl -in $QUERY -out $QUERY
        blastn -query $QUERY -db $DATABASE -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" -out $OUTPUT -num_threads $THREADS


else
    echo "Program input incorrectly formatted. Please input blastn/BLASTN, tblastn/TBLASTN or blastx/BLASTX"
fi
