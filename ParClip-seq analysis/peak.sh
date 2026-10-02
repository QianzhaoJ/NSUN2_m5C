#!/bin/bash
for sample in {RIGI,RIGI_uni}
do
echo $sample
ini=/share/home/jiqianzhao/05_Results/02_NSUN2/parclip/scripts/PARalyzer_v1_5/ini/${sample}.ini
out=/share/home/jiqianzhao/05_Results/02_NSUN2/parclip/RIG1/05_peak
mkdir -p $out
./PARalyzer 5G $ini 2>>peak.log
done
