#!/bin/bash
#This is a pipeline for RNA-seq analysis
#Requirement: conda activate daily
#Current version 1: J.Q.Z, 2021-10-27

array=( "$@" )

#Usage
if [ ! -n "$1" ]
then
  echo "********************************************************************************************"
  echo "*                     PipeRiboseq: pipeline for RNA-seq analysis.                         *"
  echo "*                            Version 2, 2021-10-27, Q.Z                                    *"
  echo "* Usage: `basename $0`                                                                     *"
  echo "*        Required (paired end data):                                                       *"
  echo "*                  -F [fastq files1]                                                       *"
  echo "*                  -f [fastq files2]                                                       *"
  echo "*                  -n [The prefix samplename]                                              *"
  echo "*                  -o [The outputpath]                                                     *"
  echo "*        Optional: -p [Number of CPUs, default=24]                                         *"
  echo "*                  -g [The file of Human Gene GTF files]                                         *"
  echo "*                  -G [The file of ERCC GTF files]                                          *"
  echo "*                  -i [The file of hisat2 human index files]                                     *"
  echo "*                  -I [The file of hisat2 ERCC index files]                                       *"
  echo "* Inputs: The raw fastq files                                                              *"
  echo "* Run: Default to run trim_galore,fastqc,hisat2,sam2bam,featureCounts & deeptools          *"
  echo "*      Figures will be generated in /plots folder, and bigWig files in /tracks folder      *"
  echo "* Outputs: All output files will be generated in the same folder as the pipeline submitted *"
  echo "* This pipeline requires 'conda activate daily'                                            *"
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
egtf="unassigned"
index="unassigned"
eindex="unassigned"

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
    echo '   Human GTF: '$gtf
 elif [[ $arg == "-G" ]]
  then
    egtf=${array[$counter+1]}
    echo '   ERCC GTF: '$egtf
 elif [[ $arg == "-i" ]]
  then
    index=${array[$counter+1]}
    echo '   Index: '$index
 elif [[ $arg == "-I" ]]
  then
    eindex=${array[$counter+1]}
    echo '   ERCC Index: '$eindex
 fi
  let counter=$counter+1
done

echo "*****************************************************************************"
echo "1. Generate the all output files"

cd $out
mkdir -p 01_fastqc 02_trim 03_mapping 04_count 05_repeat 06_FPKM

trim=$out/02_trim
map=$out/03_mapping
count=$out/04_count
repeat=$out/05_repeat
fpkm=$out/06_FPKM

for path in `ls $out | grep 0`
do
mkdir -p $path/logs
done

echo "*****************************************************************************"
echo "2. Run trim_galore"

trim_log=$trim/logs/${sample}.log
trim_out=$trim/$sample
mkdir -p $trim_out
#trim_galore --fastqc --fastqc_args "--threads 24" --cores 12 --path_to_cutadapt ~/anaconda3/envs/daily/bin/cutadapt --stringency 3 --paired --output_dir $trim_out $fq1 $fq2 2>$trim_log

echo "*****************************************************************************"
echo "3. Mapping using hisat2 Human"

map_log=$map/logs/${sample}.hisat2.log
map_out=$map/$sample
mkdir -p $map_out
# input ###
clean1=$trim_out/${sample}*1.fq.gz
clean2=$trim_out/${sample}*2.fq.gz
# out ###
sam=$map_out/${sample}.sam
bam=$map_out/${sample}.bam
sortCoord=$map_out/${sample}.sortCoord.bam

hisat2 --rna-strandness RF -p $CPU -x $index -1 $clean1 -2 $clean2 --dta -S $map_out/${sample}.sam 2>$map_log
samtools view -bS -@ $CPU -q 10 $map_out/${sample}.sam > $bam &&
rm $map_out/${sample}.sam
samtools sort -@ $CPU $bam -o $sortCoord
samtools index $sortCoord

echo "*****************************************************************************"
echo "3.2. Mapping using hisat2 ERCC"

emap_log=$map/logs/${sample}.hisat2.ERCC.log
esam=$map_out/${sample}.ERCC.sam
ebam=$map_out/${sample}.ERCC.bam
esortCoord=$map_out/${sample}.sortCoord.ERCC.bam

hisat2 --rna-strandness RF -p $CPU -x $eindex -1 $clean1 -2 $clean2 --dta -S $esam 2>$emap_log
samtools view -bS -@ $CPU -q 10 $esam > $ebam &&
rm $esam
samtools sort -@ $CPU $ebam -o $esortCoord
samtools index $esortCoord


echo "*****************************************************************************"
echo "4. assigning sequence reads to genomic features"

count_out=$count/$sample
mkdir -p $count_out

featureCounts -o $count_out/${sample}.hisat2.txt -Q 10 -p -T $CPU -F GTF -a $gtf -t exon -g gene_id -s 2 -B -C --minOverlap 10 $sortCoord
awk '{print $1,$7}' $count_out/${sample}.hisat2.txt | grep -v "#" >$count_out/${sample}.hisat2.final.txt

featureCounts -o $count_out/${sample}.hisat2.ERCC.txt -Q 10 -p -T $CPU -F GTF -a $egtf -t exon -g gene_id -s 2 -B -C --minOverlap 10 $esortCoord
awk '{print $1,$7}' $count_out/${sample}.hisat2.ERCC.txt | grep -v "#" >$count_out/${sample}.hisat2.ERCC.final.txt
