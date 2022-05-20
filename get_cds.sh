#!/bin/bash

 zcat ./data/gencode.v19.annotation.gff3.gz | grep -w "$1" |\
  agrep -w "CCDS;CDS" |  grep -v "transcript_type=retained_intron" |\
   grep -P "appris_principal|appris_candidate_longest" |\
    awk -v FS="\t" -v OFS="\t" '{split($9, a, ";"); print $1, $2, a[1]":"a[7]":"$7, $4-1, $5+1, $6, $7, $8, $9}'

