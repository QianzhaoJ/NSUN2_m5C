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
  echo "*                  -g [The file of Gene GTF files]                                         *"
  echo "*                  -G [The file of REs GTF files]                                          *"
  echo "*                  -i [The file of hisat2 index files]                                     *"
  echo "*                  -I [The file of STAR index files]                                       *"
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
REgtf="unassigned"
index="unassigned"
STARindex="unassigned"

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
 elif [[ $arg == "-G" ]]
  then
    REgtf=${array[$counter+1]}
    echo '   REGTF: '$REgtf
 elif [[ $arg == "-i" ]]
  then
    index=${array[$counter+1]}
    echo '   Index: '$index
 elif [[ $arg == "-I" ]]
  then
    STARindex=${array[$counter+1]}
    echo '   STARIndex: '$STARindex
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
#mkdir -p  $path/scripts
mkdir -p $path/logs
done

echo "*****************************************************************************"
echo "2. Run trim_galore"

trim_log=$trim/logs/${sample}.log
trim_out=$trim/$sample
mkdir -p $trim_out
trim_galore --fastqc --path_to_cutadapt ~/anaconda3/envs/daily/bin/cutadapt --stringency 3 --paired --output_dir $trim_out $fq1 $fq2 2>$trim_log

echo "*****************************************************************************"
echo "3. Mapping using hisat2"

map_log=$map/logs/${sample}.log
map_out=$map/$sample
mkdir -p $map_out
# input ###
clean1=$trim_out/${sample}*1.fq.gz
clean2=$trim_out/${sample}*2.fq.gz
# out ###
bam=$map_out/${sample}.bam
sortCoord=$map_out/${sample}.sortCoord.bam
sortName=$map_out/${sample}.sortName.bam
bw=$map_out/${sample}.bw
STARout=$map_out/STAR

hisat2 -p $CPU -x $index -1 $clean1 -2 $clean2 --dta -S $map_out/${sample}.sam 2>$map_log
samtools view -bS -@ $CPU -q 10 $map_out/${sample}.sam > $bam &&
rm $map_out/${sample}.sam
samtools sort -@ $CPU $bam -o $sortCoord
samtools sort -@ $CPU -n $bam -o $sortName
samtools index $sortCoord

bamCoverage -p $CPU -b $sortCoord -o $bw --normalizeUsing RPKM --binSize 10 --ignoreForNormalization chrY

echo "*****************************************************************************"
echo "3.2. Mapping using STAR"

mkdir -p $STARout

#STAR --genomeDir $STARindex --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --runThreadN $CPU --readFilesIn $clean1 $clean2 --outSAMunmapped None --outFileNamePrefix $STARout/${sample}.STAR. --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --outFilterMismatchNoverLmax 0.04 --outFilterType BySJout 

#STAR --runThreadN 24 --genomeDir $STARindex --outSAMtype BAM SortedByCoordinate --runMode alignReads --outFilterMultimapNmax 1000 --outSAMmultNmax 1 --outFilterMismatchNmax 3 --outMultimapperOrder Random --winAnchorMultimapNmax 1000 --alignEndsType EndToEnd --alignIntronMax 1--alignMatesGapMax 350 --readFilesIn $clean1 $clean2 --readFilesCommand zcat --outFileNamePrefix $STARout/${sample}.STAR.
echo "*****************************************************************************"
echo "4. assigning sequence reads to genomic features"

count_log=$count/logs/${sample}.log
count_out=$count/$sample
mkdir -p $count_out

featureCounts -T $CPU -Q 10 -p -t exon -g gene_id -a $gtf -o $count_out/${sample}.rna.txt $sortCoord
awk '{print $1,$7}' $count_out/${sample}.rna.txt | grep -v "#" >$count_out/${sample}_final.txt

echo "*****************************************************************************"
echo "5. assigning sequence reads to repeat elements"

STARbam=$STARout/${sample}.STAR.Aligned.sortedByCoord.out.bam
repeat_log=$repeat/logs/${sample}.log
repeat_out=$repeat/$sample
mkdir -p $repeat_out

#featureCounts -M -T $CPU -p -t exon -g gene_id -a $REgtf -o $repeat_out/${sample}.RE.txt $STARbam

echo "*****************************************************************************"
echo "6. Quantitation of RNA-Seq data using stringtie"

stringtie=/share/home/jiqianzhao/anaconda3/envs/FPKM/bin/stringtie
fpkm_log=$fpkm/logs/${sample}.log
fpkm_out=$fpkm/$sample
fpkm_all=$fpkm_out/${sample}.gene_abund.tab
fpkm_gtf=$fpkm_out/${sample}.gtf

$stringtie -A $fpkm_all -e -B -p $CPU -G $gtf -o $fpkm_gtf $sortCoord 2>$fpkm_log
