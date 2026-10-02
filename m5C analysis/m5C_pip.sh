#!/bin/bash
#This is a pipeline for m5C analysis
#Requirement: conda activate m5C
#Current version 1: J.Q.Z, 2022-05-06

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                     PipeRiboseq: pipeline for m5C-seq analysis.                         *"
  echo "*                            Version 2, 2022-05-06, Q.Z                                    *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -F [fastq files1]                                                       *"
  echo "*                  -f [fastq files2]                                                       *"
  echo "*                  -n [The prefix samplename]                                              *"
  echo "*                  -o [The outputpath]                                                     *"
  echo "*        Optional: -p [Number of CPUs, default=24]                                         *"
  echo "*                  -g [The file of Gene GTF files]                                         *"
  echo "*                  -i [The file of meRanT index files]                                     *"
  echo "* Inputs: The raw fastq files                                                              *"
  echo "* Run: Default to run trim_galore,fastqc,meRanGh,meRanCall                                 *"
  echo "* Outputs: All output files will be generated in outputpath *"
  echo "* This pipeline requires 'conda activate m5C'                                              *"
  echo "********************************************************************************************"
  exit 1
fi

#Get parameters
fq1="unassigned"
fq2="unassigned"
sample="unassigned"
out="unassigned"
CPU=24
gtf="unassigned"
index="unassigned"

for arg in "$@"
do
 if [[ $arg == "-F" ]]
  then
    fq1=${array[$counter+1]}
    echo '   raw fastq1: '$fq1
 elif [[ $arg == "-f" ]]
  then
    fq2=${array[$counter+1]}
    echo '   raw fastq2: '$fq2
 elif [[ $arg == "-n" ]]
  then
    sample=${array[$counter+1]}
    echo '   Sample name: '$sample
 elif [[ $arg == "-o" ]]
  then
    out=${array[$counter+1]}
    echo '   The output of all files: '$out
 elif [[ $arg == "-p" ]]
  then
    CPU=${array[$counter+1]}
    echo '   CPU: '$CPU
 elif [[ $arg == "-g" ]]
  then
    gtf=${array[$counter+1]}
    echo '   GTF: '$gtf
 elif [[ $arg == "-i" ]]
  then
    index=${array[$counter+1]}
    echo '   Index: '$index
 fi
  let counter=$counter+1
done

echo "*****************************************************************************"
echo "1. Generate the all output files"

cd $out
mkdir -p 01_fastqc 02_trim 03_mapping 04_DHFR 05_methyC

trim=$out/02_trim
map=$out/03_mapping
dhfr=$out/04_DHFR
call=$out/05_methyC

for path in `ls $out | grep 0`
do
mkdir -p $path/logs
done

echo "*****************************************************************************"
echo "2. Run trim_galore"

trim_log=$trim/logs/${sample}.log
trim_out=$trim/$sample
mkdir -p $trim_out

/share/home/jiqianzhao/anaconda3/envs/daily/bin/cutadapt -a AGATCGGAAGAGCACACGTCTG -A AGATCGGAAGAGCGTCGTGTAG -o ${trim_out}/${sample}_1_val_1.fq.gz -p ${trim_out}/${sample}_2_val_2.fq.gz $fq1 $fq2 2>$trim_log

trimmomatic PE -phred33 -threads 24 ${trim_out}/${sample}_1_val_1.fq.gz ${trim_out}/${sample}_2_val_2.fq.gz ${trim_out}/${sample}_1.trim.fq.gz ${trim_out}/${sample}_1.unpaired.fq.gz ${trim_out}/${sample}_2.trim.fq.gz ${trim_out}/${sample}_2.unpaired.fq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:30 MINLEN:35 2>>$trim_log

/share/home/jiqianzhao/anaconda3/envs/daily/bin/fastqc -t 12 -o ${trim_out} ${trim_out}/${sample}_1.trim.fq.gz ${trim_out}/${sample}_2.trim.fq.gz

echo "*****************************************************************************"
echo "3. Mapping using meRanG"

map_log=$map/logs/${sample}.log
map_out=$map/$sample
mkdir -p $map_out
# input ###
clean1=${trim_out}/${sample}_1.trim.fq.gz
clean2=${trim_out}/${sample}_2.trim.fq.gz
# out ###
sam=$map_out/${sample}.sam

/share/home/jiqianzhao/04_Softwares/meRanTK-1.2.1b/meRanGh align -o $map_out -f $clean2 -r $clean1 -t 12 -S $sam -ds -id $index -GTF $gtf -fmo -mmr 0.01 2>$map_log

echo "*****************************************************************************"
echo "4. Mapping using meRanT"

dhfr_log=$dhfr/logs/${sample}.log
dhfr_out=$dhfr/$sample
mkdir -p $dhfr_out
# input ###
clean1=$trim_out/${sample}*1.fq.gz
clean2=$trim_out/${sample}*2.fq.gz
# out ###
DHFRsam=$dhfr_out/${sample}.DHFR.meRanTK.sam

/share/home/jiqianzhao/04_Softwares/meRanTK-1.2.1b/meRanT align -o $dhfr_out -f $clean2 -r $clean1 -t 12 -k 10 -S $DHFRsam -un -ds -i2g /share/home/jiqianzhao/03_Database/m5C/DHFR/DHFR.map -x /share/home/jiqianzhao/03_Database/m5C/DHFR/DHFR/DHFR940-1503_C2T -mbp -fmo -mmr 0.01 2>$dhfr_log

/share/home/jiqianzhao/04_Softwares/meRanTK-1.2.1b/meRanCall -p 12 -f /share/home/jiqianzhao/03_Database/m5C/DHFR/DHFR940-1503.fa -bam $dhfr_out/${sample}.DHFR.meRanTK_sorted.bam -o $dhfr_out/${sample}.meRanCall.result -tref -mBQ 30 -mr 0 >>$dhfr_log

sed '1d' $dhfr_out/*.meRanCall.result | awk '{C=C+$5; m5C+=$6}END{print (C-m5C)/C}' > $dhfr_out/CT_conversion.txt

echo "*****************************************************************************"
echo "5. meRanCall"

call_log=$call/logs/${sample}.log
call_out=$call/$sample
mkdir -p $call_out

/share/home/jiqianzhao/04_Softwares/meRanTK-1.2.1b/meRanCall -p 12 -f /share/home/jiqianzhao/03_Database/m5C/Homo_sapiens.GRCh37.dna_sm.primary_assembly.fa  -bam $map_out/${sample}_sorted.bam  -gref -o $call_out/${sample}.fdr.txt -mBQ 20 -mr 0 -fdr 0.05 2>$call_log
cat $call_out/${sample}.fdr.txt | awk -v OFS="\t" '{if($12=="M")print $1,$2,$2,$3,$6,$5,$7,$12,$15,$19}' > $call_out/${sample}.methyC.fdr.txt
cat $call_out/${sample}.methyC.fdr.txt | awk -v OFS="\t" '$7>0 && $6>0 && $5>0' > $call_out/${sample}.methyC_qua0.txt 
