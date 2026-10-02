#!/bin/bash

tar=/share/home/jiqianzhao/05_Results/02_NSUN2/parclip
raw=/share/home/jiqianzhao/02_Data/NSUN2/parclip
trim=$tar/02_trim
map=$tar/04_True

mkdir -p $map/scripts

for sample in `ls $raw`
do
echo -e '#!/bin/bash' > $map/scripts/${sample}_map.sh
echo '#SBATCH -N 1' >> $map/scripts/${sample}_map.sh
echo "#SBATCH -n 1" >> $map/scripts/${sample}_map.sh
echo "#SBATCH -c 24" >> $map/scripts/${sample}_map.sh
echo "#SBATCH -p compute" >> $map/scripts/${sample}_map.sh
echo "#SBATCH --job-name=Map_$sample" >> $map/scripts/${sample}_map.sh
echo '#SBATCH --export=ALL' >> $map/scripts/${sample}_map.sh
echo -e tar=$tar\\n >> $map/scripts/${sample}_map.sh
echo -e sample=${sample}\\ntar=$tar\\nraw=$raw'

echo Mapping: $sample is Starting
gtf=/share/home/jiqianzhao/03_Database/m5C/Homo_sapiens.GRCh37.87.gtf
trim=$tar/02_trim
map=$tar/04_True
log=$map/logs

results=$map/$sample
mkdir -p $results

clean=$trim/$sample/${sample}_R1.trimmed.uni.fq

genome=/share/home/jiqianzhao/03_Database/m5C/m5C_bowtie/hg19

bowtie $genome -p 24 -v 2 -m 10 --best --strata -S $clean $results/${sample}.bowtie_true.sam
htseq-count -m union -s reverse $results/${sample}.bowtie_true.sam $gtf > $results/${sample}.count.txt
###############
/share/home/jiqianzhao/anaconda3/envs/daily/bin/samtools view -@ 12 -b -q 20 $results/${sample}.bowtie_true.sam > $results/${sample}.bowtie_true.bam
/share/home/jiqianzhao/anaconda3/envs/daily/bin/samtools sort -@ 12 $results/${sample}.bowtie_true.bam > $results/${sample}.bowtie_true.sorted.bam
/share/home/jiqianzhao/anaconda3/envs/daily/bin/samtools index $results/${sample}.bowtie_true.sorted.bam

echo Mapping has been Done' >> $map/scripts/${sample}_map.sh
sbatch $map/scripts/${sample}_map.sh
done
