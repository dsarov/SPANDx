#!/bin/bash


##This script is run as part of the SPANDx pipeline as of version 3.1

PBS_O_WORKDIR=$1
ref=$2
SAMTOOLS=$3

##input sample names

cd $PBS_O_WORKDIR

sequences_tmp=(`find $PBS_O_WORKDIR/*_1_sequence.fastq.gz -printf "%f "`)
sequences=("${sequences_tmp[@]/_1_sequence.fastq.gz/}")
n=${#sequences[@]}

TAB="$(printf '\t')"

cat << _EOF_ > SNP_summary_header

Reference=$ref
Total number of genomes analysed=$n
Genome Name${TAB}Number of SNPs passing filters${TAB}Number of SNPs Failing filters${TAB}Number of indels passing filters${TAB}Number of indels failing filters${TAB}Average coverage${TAB}Mapped reads
_EOF_



for (( i=0; i<n; i++ )); do
  
  SNP_PASS_count=$(cat Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf | grep -v '#' | wc -l)
  SNP_FAIL_count=$(cat Outputs/SNPs_indels_FAIL/${sequences[$i]}.snps.FAIL.vcf | grep -v '#' | wc -l)
  indel_PASS_count=$(cat Outputs/SNPs_indels_PASS/${sequences[$i]}.indels.PASS.vcf | grep -v '#' | wc -l)
  indel_FAIL_count=$(cat Outputs/SNPs_indels_FAIL/${sequences[$i]}.indels.FAIL.vcf | grep -v '#' | wc -l)
  Avg_cov=$(cat ${sequences[$i]}/unique/${sequences[$i]}.sample_summary | awk '{print $3}' | head -n2 | tail -n1)
  Mapped_reads=$($SAMTOOLS idxstats ${sequences[$i]}/unique/${sequences[$i]}.realigned.bam | awk '{ print $3 }' | head -n1)
  
  echo -e "${sequences[$i]}\t$SNP_PASS_count\t$SNP_FAIL_count\t$indel_PASS_count\t$indel_FAIL_count\t$Avg_cov\t$Mapped_reads" > ${sequences[$i]}.summary
  
done  
  
cat SNP_summary_header ${sequences[@]/%/.summary} > Outputs/Single_sample_summary.txt

rm SNP_summary_header
rm ${sequences[@]/%/.summary} 

exit 0
