# Usage

```bash
pip=/share/home/jiqianzhao/01_scripts/m5C/m5C_pip.sh
fq1=\(raw/\)sample/${sample}_*1.fq.gz
fq2=\(raw/\)sample/${sample}_*2.fq.gz

$pip -F $fq1 -f $fq2 -n $sample -o $tar -p 12 -g $gtf -i $index
```
