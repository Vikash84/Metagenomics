#!/bin/bash

set -e
set -u
set -o pipefail

samples=${1?missing sample file}
threads=8
OUT_DIR="analysis/2-read-align"

echo -e "Experiment \t Raw Reads \t Uniquely mapped Reads \t ratio"
exec 0< $samples

while read sample;
do
     total=$( samtools view -@ ${threads} -c ${ALIGN_DIR}/${sample}.sort.bam )
     unique=$( samtools view -@ ${threads} -c ${ALIGN_DIR}/${sample}.flt.bam )
     ratio=$( echo "scale=2; 100 * $unique / $total " | bc )
     echo -e "$sample \t $total \t $unique \t $ratio %"
done




Experiment	Raw Reads	Uniquely mapped Reads	ratio
HKE293-D210N-V5ChIP-Rep1	22405416	6443979	28.76%
HKE293-D210N-Input-Rep1	60302237	25673307	42.57%
HKE293-D210N-V5ChIP-Rep2	17763614	11778533	66.30%
HKE293-D210N-Input-Rep2	11131443	8553097	76.83%
HKE293-D210N-V5ChIP-Rep3	8799855	5640375	64.09%
HKE293-D210N-Input-Rep3	4529910	3209275	70.84%
HKE293-WKKD-V5ChIP	12734577	8612940	67.63%
HKE293-WKKD-Input	8830478	6643507	75.23%
HKE293-delta-HC-V5ChIP	25174573	9252009	36.75%
