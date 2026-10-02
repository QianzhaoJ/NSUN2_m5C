#!/bin/bash

tar=/share/home/jiqianzhao/05_Results/02_NSUN2/parclip
raw=/share/home/jiqianzhao/02_Data/NSUN2/parclip
trim=$tar/02_trim

for sample in `ls $raw`
do
echo -e '#!/bin/bash' > $trim/scripts/${sample}_trim.sh
echo '#SBATCH -N 1' >> $trim/scripts/${sample}_trim.sh
echo "#SBATCH -n 1" >> $trim/scripts/${sample}_trim.sh
echo "#SBATCH -c 12" >> $trim/scripts/${sample}_trim.sh
echo "#SBATCH -p compute" >> $trim/scripts/${sample}_trim.sh
echo "#SBATCH --job-name=Trim_$sample" >> $trim/scripts/${sample}_trim.sh
echo '#SBATCH --export=ALL' >> $trim/scripts/${sample}_trim.sh
echo -e tar=$tar\\n >> $trim/scripts/${sample}_trim.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Trimming: $sample is Starting

trim=$tar/02_trim
log=$trim/logs

results=$trim/$sample
mkdir -p $results

fq=$raw/$sample/${sample}_R1.fq.gz
fq1=$raw/$sample/${sample}_R1.fq
step1=$results/${sample}_R1.step1.fq
step2=$results/${sample}_R1.step2.fq
step3=$results/${sample}_R1.step3.fq
final=$results/${sample}_R1.trimmed.fq
uni=$results/${sample}_R1.trimmed.uni.fq

gunzip $fq
cutadapt -a AGATCGGAAGAGCACACGTCTG -o $step1 $fq1
cutadapt -a AAAAAAAA -o $step2 $step1
cutadapt -u 4 -o $step3 $step2
java -Xmx40g -jar /share/home/jiqianzhao/04_Softwares/Trimmomatic-0.39/trimmomatic-0.39.jar SE -phred33 $step3 $final LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:15
perl /share/home/jiqianzhao/01_scripts/RICseq/RICpipe-master/step0.remove_PCR_duplicates/scripts/remove_duplicated_reads_SE.pl $final $uni
fastqc $step1 -t 12 -o $results
fastqc $step3 -t 12 -o $results
fastqc $final -t 12 -o $results
fastqc $uni -t 12 -o $results
echo Trimming has been Done' >> $trim/scripts/${sample}_trim.sh
sbatch $trim/scripts/${sample}_trim.sh
done

