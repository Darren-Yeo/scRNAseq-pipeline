

#!/bin/bash


files2=("CRR1984896" "CRR1984894" "CRR1984897")


for f2 in "${files2[@]}"
do
  wget -c ftp://download.big.ac.cn/gsa6/CRA027693/${f2}/${f2}_r1.fastq.gz
  wget -c ftp://download.big.ac.cn/gsa6/CRA027693/${f2}/${f2}_r2.fastq.gz
done
