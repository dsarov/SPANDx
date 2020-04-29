#!/bin/bash

#########################################################################
# The following script will convert SNP and/or indel calls into a nexus file and create a preliminary tree with Fastree
#
#
# Written by D. Sarovich
# dsarovich@usc.edu.au
#
#########################################################################

#set variant genome
#variant_genome_path=$1
baseDir=$1

echo "Creating VCF tables"
gatk VariantsToTable -V out.filtered.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.vcf.table
gatk VariantsToTable -V out.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.vcf.table.all


###########################################################################
# Creates the SNP matrix for PAUP
###########################################################################

echo "Creating SNP matrix"
awk '$5 ~/SNP/' out.vcf.table | awk '$4 !~/\.,\*/' | grep -v '\./\.' > out.vcf.table.snps.clean
# replace A/A, C/C, G/G, T/T genotypes with single nucleotides A, G, C, T etc etc 
sed -i 's#A/A\|A|A#A#g' out.vcf.table.snps.clean
sed -i 's#G/G\|G|G#G#g' out.vcf.table.snps.clean
sed -i 's#C/C\|C|C#C#g' out.vcf.table.snps.clean
sed -i 's#T/T\|T|T#T#g' out.vcf.table.snps.clean
grep -v '|' out.vcf.table.snps.clean | grep -v '/' > vcf.table.tmp #remove mixed genotypes
mv vcf.table.tmp out.vcf.table.snps.clean
taxa=$(head -n1 out.vcf.table | cut -f3,6- | sed 's/.GT//g')
ntaxa=$(awk '{print NF-4; exit }' out.vcf.table)
nchar=$(cat out.vcf.table.snps.clean | wc -l)
awk '{print $1,$2}' out.vcf.table.snps.clean | sed 's/ /_/g' > snp.location
cut -f3,6- out.vcf.table.snps.clean > grid.nucleotide
grid=$(paste snp.location grid.nucleotide)
echo -e "\n#nexus\nbegin data;\ndimensions ntax=$ntaxa nchar=$nchar;\nformat symbols=\"AGCT\" gap=. transpose;\ntaxlabels $taxa;\nmatrix\n$grid\n;\nend;" > Ortho_SNP_matrix.nex

###########################################################################
# Creates the SNP matrix for FastTree2
###########################################################################

awk 'BEGIN {FS=" "; OFS=""}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' grid.nucleotide > grid.nucleotide.fasttree
head -n1 out.vcf.table | cut -f3,6- | sed 's/.GT//g' > taxa.tmp
awk 'BEGIN {FS=" "; OFS=""}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' taxa.tmp > taxa.fasttree
sed -i 's/^/>/' taxa.fasttree
paste -d '\n' taxa.fasttree grid.nucleotide.fasttree >Ortho_SNP_matrix_FastTree2.nex
fasttree -log ML_log.txt -nt Ortho_SNP_matrix_FastTree2.nex > ML_phylogeny.tre
fasttree -log MP_log.txt -noml -nt Ortho_SNP_matrix_FastTree2.nex > MP_phylogeny.tre


	
echo "SPANDx has finished"

exit 0
