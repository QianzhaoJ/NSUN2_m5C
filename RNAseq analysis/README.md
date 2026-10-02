# Usage

```bash
pip=RNASeq_pip.sh
fq1=$raw/$sample/${sample}_*1.fq.gz
fq2=$raw/$sample/${sample}_*2.fq.gz

$pip -F $fq1 -f $fq2 -n $sample -o $tar -p 24 -g $gtf -G $REgtf -i $index -I $STARindex
```
