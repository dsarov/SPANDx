#!/bin/bash

#########################################################################
# The following script will convert SNP and/or indel calls into a nexus file and create a preliminary tree with Paup
# Additionally, the SNPs and/or indels are annotated with SnpEff
#
# Written by D. Sarovich
# dsarovich@usc.edu.au
#
#########################################################################

#set variant genome
variant_genome_path=$1
baseDir=$2

echo "Creating VCF tables"
gatk VariantsToTable -V out.filtered.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.vcf.table
gatk VariantsToTable -V out.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -O out.vcf.table.all


###########################################################################
# Creates the SNP matrix for PAUP
###########################################################################

echo "Creating SNP matrix"
grep 'SNP' out.vcf.table | grep -v ',' | grep -v '*' | grep -v '\.' > out.vcf.table.snps.clean
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

##run paup to create tree

#if [ "$ntaxa" -gt 4 ]; then
#  echo "Running PAUP"
#PAUP block to be inserted into nexus file
#cat <<_EOF_ > tmpnex
#begin paup;
#Set AllowPunct=Yes;
#lset nthreads=2;
#hsearch;
#savetrees from=1 to=1 brlens=yes;
#_EOF_
#  cat Ortho_SNP_matrix.nex tmpnex > run.nex
#  paup -n run.nex >& paup_log.txt
#  mv run.tre MP_phylogeny.tre
#else
#  echo "Fewer than 4 taxa found. Skipping creation of MP tree"
#fi

###############################################
##
## These steps will take the merged SNP outputs, annotate them and reformat the data into a tab delimited txt file for importation into excel
## output from these steps should be a tab delimited file
##
## Annotation of SNPs using snpEff
##
##################################################

#clean-up the out.vcf.table.all because GATK outputs A/A
sed -i 's#|#/#g' out.vcf.table.all
awk ' { for (i=6; i<=NF; i++) {
        if ($i == "A/A") $i="A"; 
        if ($i == "G/G") $i="G"; 
        if ($i == "C/C") $i="C"; 
        if ($i == "T/T") $i="T"; 
        if ($i == "*/*") $i="*"; 
        if ($i == "./.") $i=".";
        }};
        {print $0} ' out.vcf.table.all > out.vcf.table.all.tmp	 
		 
awk ' { for (i=6; i<=NF; i++) {
        if ($i ~ /\//) { 
          split($i, a, "/");
        if (a[1] == a[2]) $i=a[1];
           }
         }
       }; 
       {print $0} ' out.vcf.table.all.tmp > out.vcf.table.all

	
snpEff eff -no-downstream -no-intergenic -ud 100 -formatEff -v -dataDir ${baseDir}/resources/snpeff ${variant_genome_path} out.vcf > out.annotated.vcf
	
#remove headers from annotated vcf and out.vcf
grep -v '#' out.annotated.vcf > out.annotated.vcf.headerless
#grep -v '#' out.vcf > out.vcf.headerless
awk '{
    if (match($0,"EFF=")){print substr($0,RSTART)}
    else
    print ""
    }' out.annotated.vcf.headerless > effects

sed -i 's/EFF=//' effects
sed -i 's/(/ /g' effects
sed -i 's/|/ /g' effects
sed -i 's/UPSTREAM MODIFIER /UPSTREAM MODIFIER - /g' effects
cut -d " " -f -8 effects > effects.mrg
sed -i 's/ /\t/g' effects.mrg
rm effects
    
tail -n+2 out.vcf.table.all > out.vcf.table.all.headerless
sed -i 's/ /\t/g' out.vcf.table.all.headerless
paste out.vcf.table.all.headerless effects.mrg > out.vcf.headerless.plus.effects
head -n1 out.vcf.table.all | sed 's/.GT//g' > header.left
echo -e "Effect\tImpact\tFunctional_Class\tCodon_change\tProtein_and_nucleotide_change\tAmino_Acid_Length\tGene_name\tBiotype" > header.right
paste header.left header.right > header
cat header out.vcf.headerless.plus.effects > All_SNPs_indels_annotated.txt
	
echo "ARDaP has finished"

exit 0
