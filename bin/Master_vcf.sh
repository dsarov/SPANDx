#!/bin/bash

ref=$1

echo "Creating master VCF file"

array=($(find *.gvcf -printf "%f "))
n="${#array[@]}"
array2=("${array[@]/#/-V }")
gatk CombineGVCFs -R ${ref}.fasta ${array2[*]} -O master.vcf
gatk GenotypeGVCFs -ploidy $n -R ${ref}.fasta -V master.vcf -O out.vcf
 
exit 0
