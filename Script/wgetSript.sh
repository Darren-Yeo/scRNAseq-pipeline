
#!/bin/bash


files1=("CRR1984893" "CRR1984898" "CRR1984895")



for f in "${files1[@]}"
do
  wget -c ftp://download.big.ac.cn/gsa6/CRA027693/${f}/${f}_r1.fastq.gz
  wget -c ftp://download.big.ac.cn/gsa6/CRA027693/${f}/${f}_r2.fastq.gz
done


### example of what it does
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984893/CRR1984893_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984893/CRR1984893_r2.fastq.gz
#
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984898/CRR1984898_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984898/CRR1984898_r2.fastq.gz
#
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984895/CRR1984895_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984895/CRR1984895_r2.fastq.gz
#
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984896/CRR1984896_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984896/CRR1984896_r2.fastq.gz
#
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984894/CRR1984894_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984894/CRR1984894_r2.fastq.gz
#
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984897/CRR1984897_r1.fastq.gz
#wget ftp://download.big.ac.cn/gsa6/CRA027693/CRR1984897/CRR1984897_r2.fastq.gz

exit