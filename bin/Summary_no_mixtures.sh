#!/bin/bash



ref=$1
baseDir=$2

##input sample names

#cd $PBS_O_WORKDIR

sequences_tmp=(`find ../../../Outputs/Variants/VCFs/*.PASS.snps.vcf -printf "%f "`)
sequences=("${sequences_tmp[@]/.PASS.snps.vcf/}")
n=${#sequences[@]}

TAB="$(printf '\t')"

cat << _EOF_ > SNP_summary_header

Reference=$ref
Total number of genomes analysed=$n
Genome Name${TAB}Number of SNPs passing filters${TAB}Number of SNPs failing filters${TAB}Number of indels passing filters${TAB}Number of indels failing filters${TAB}Number of mixed variants${TAB}Average coverage${TAB}Mapped reads
_EOF_


#Mix summary
read -r HEADER < All_SNPs_indels_annotated.txt
NUM_COLS=$(awk '{print NF}' <<<"$HEADER" FS="\t")

# Print the header for the output
echo -e "Sample Name\tCount" > mixture_sum.txt

# Loop through each column by its index
for ((i=6; i<=NUM_COLS-8; i++))
do
    # Extract the column header
    HEADER_NAME=$(cut -f $i <<<"$HEADER")
    
    # Count the number of entries with "/" in the current column
    COUNT=$(awk -v col="$i" -F'\t' 'NR>1 && $col ~ /\// {count++} END {print count+0}' All_SNPs_indels_annotated.txt)
    
    # Print the result
    echo -e "${HEADER_NAME}\t${COUNT}" >> mixture_sum.txt
done



for (( i=0; i<n; i++ )); do
  
  snp_PASS_count=$(cat ../../../Outputs/Variants/VCFs/${sequences[$i]}.PASS.snps.vcf | grep -v '#' | wc -l)
  snp_FAIL_count=$(cat ../../../Outputs/Variants/VCFs/${sequences[$i]}.FAIL.snps.vcf | grep -v '#' | wc -l)
  indel_PASS_count=$(cat ../../../Outputs/Variants/VCFs/${sequences[$i]}.PASS.indels.vcf | grep -v '#' | wc -l)
  indel_FAIL_count=$(cat ../../../Outputs/Variants/VCFs/${sequences[$i]}.FAIL.indels.vcf | grep -v '#' | wc -l)
  Mapped_reads=$(samtools idxstats ../../../Outputs/bams/${sequences[$i]}.dedup.bam |  awk '{ total += $3 } END {print total }')
  Avg_cov=$(cat "${sequences[$i]}".depth.txt)
  Mix_count=$(grep -w "${sequences[$i]}" mixture_sum.txt | awk '{print $2}')
  #Mixed_variants=$(All_SNPs_indels_annotated.txt
  
  
  echo -e "${sequences[$i]}\t$snp_PASS_count\t$snp_FAIL_count\t$indel_PASS_count\t$indel_FAIL_count\t$Mix_count\t$Avg_cov\t$Mapped_reads" > ${sequences[$i]}.summary
 
  
done
  
cat SNP_summary_header ${sequences[@]/%/.summary} > QC_metrics_summary.tsv

#rm SNP_summary_header
#rm ${sequences[@]/%/.summary} 


#exit 0
