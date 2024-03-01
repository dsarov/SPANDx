#!/bin/bash



ref=$1
baseDir=$2

##input sample names

#cd $PBS_O_WORKDIR

sequences_tmp=(`find $baseDir/Outputs/Variants/VCFs/*.PASS.snps.indels.mixed.vcf -printf "%f "`)
sequences=("${sequences_tmp[@]/.PASS.snps.indels.mixed.vcf/}")
n=${#sequences[@]}

TAB="$(printf '\t')"

cat << _EOF_ > SNP_summary_header

Reference=$ref
Total number of genomes analysed=$n
Genome Name${TAB}Number of SNPs passing filters${TAB}Number of SNPs Failing filters${TAB}Number of indels passing filters${TAB}Number of indels failing filters${TAB}Average coverage${TAB}Mapped reads
_EOF_



for (( i=0; i<n; i++ )); do
  
  variant_PASS_count=$(cat $baseDir/Outputs/Variants/VCFs/${sequences[$i]}.PASS.indels.snps.mixed.vcf | grep -v '#' | wc -l)
  variant_FAIL_count=$(cat $baseDir/Outputs/SNPs_indels_FAIL/${sequences[$i]}.FAIL.snps.indels.mixed.vcf | grep -v '#' | wc -l)
  #Avg_cov=$(cat ${sequences[$i]}/unique/${sequences[$i]}.sample_summary | awk '{print $3}' | head -n2 | tail -n1)
  Mapped_reads=$(samtools idxstats $baseDir/Outputs/bams/${sequences[$i]}.dedup.bam |  awk '{ total += $3 } END {print total }')
  
  echo -e "${sequences[$i]}\t$SNP_PASS_count\t$SNP_FAIL_count\t$indel_PASS_count\t$indel_FAIL_count\t$Avg_cov\t$Mapped_reads" > ${sequences[$i]}.summary
  
done
  
cat SNP_summary_header ${sequences[@]/%/.summary} > Outputs/Single_sample_summary.txt

rm SNP_summary_header
rm ${sequences[@]/%/.summary} 



QC_metrics_summary.tsv


exit 0
