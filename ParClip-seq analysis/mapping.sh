#!/bin/bash

tar=/share/home/jiqianzhao/05_Results/02_NSUN2/parclip
raw=/share/home/jiqianzhao/02_Data/NSUN2/parclip
trim=$tar/02_trim
map=$tar/03_mapping

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

trim=$tar/02_trim
map=$tar/03_mapping
log=$map/logs

results=$map/$sample
mkdir -p $results

clean=$trim/$sample/${sample}_R1.trimmed.uni.fq

genome=/share/home/jiqianzhao/03_Database/m5C/m5C_bowtie/hg19
bowtie $genome -p 24 -v 2 -m 10 --best --strata $clean $results/${sample}.bowtie2.sam 
grep -v "chrMT" $results/${sample}.bowtie2.sam | awk '"'"'length($6) < 100 {print $0}'"'"' > $results/${sample}.bowtie_rmMT.sam

###############
echo '"'"'mapped_reads:'"'"' > $results/ct_conversion.info
wc $results/${sample}.bowtie_rmMT.sam >> $results/ct_conversion.info
echo '"'"''"'"'  >> $results/ct_conversion.info
echo '"'"'mutated reads:'"'"' >> $results/ct_conversion.info
cat $results/${sample}.bowtie_rmMT.sam |cut -f8|awk '"'"'$1'"'"'|wc >> $results/ct_conversion.info
echo '"'"''"'"'  >> $results/ct_conversion.info
echo '"'"'mutated nucleotide:'"'"' >> $results/ct_conversion.info
cat $results/${sample}.bowtie_rmMT.sam |cut -f8|awk '"'"'$1'"'"'|awk '"'"'{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'"'"'|awk -v FS=":" '"'"'{print $2}'"'"'|sort|uniq -c|sed -e '"'"'s/^ *//'"'"'|sed -e '"'"'s/ /\t/'"'"'  >> $results/ct_conversion.info
echo ''  >> $results/ct_conversion.info
echo '"'"'conversion ratio:'"'"' >> $results/ct_conversion.info
cat $results/${sample}.bowtie_rmMT.sam |cut -f8|awk '"'"'$1'"'"'|awk '"'"'{split($1,a,",");for(i=1;i<=length(a);i++){print a[i]}}'"'"'|awk -v FS=":" '"'"'{print $2}'"'"'|sort|uniq -c|sed -e '"'"'s/^ *//'"'"'|sed -e '"'"'s/ /\t/'"'"'|awk -v OFS="\t" '"'"'{num+=$1;tmp[$2]=$1}END{for(i in tmp){print i,tmp[i],num,tmp[i]/num}}'"'"'|sort  >> $results/ct_conversion.info


echo Mapping has been Done' >> $map/scripts/${sample}_map.sh
sbatch $map/scripts/${sample}_map.sh
done

