#!/bin/bash

#########################################################################
# The following script will convert SNP and/or indel calls into a nexus file and create a preliminary tree with FastTree
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
echo "Removing mixed SNPs, "
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
# Creates the indel matrix for PAUP
###########################################################################

echo "Creating indel matrix"
echo "Removing mixed indels"
awk '$5 ~/INDEL/' out.vcf.table | awk '$4 !~/\.,\*/' | grep -v '\./\.' | grep -v ',\*' > out.vcf.table.indels.clean
awk ' { for (i=6; i<=NF; i++) {
        if ($i ~ /\//) { 
          split($i, a, "/");
        if (a[1] == a[2]) $i=a[1];
           }
         }
       }; 
       {print $0} ' out.vcf.table.indels.clean > out.vcf.table.indels.clean.tmp	
mv out.vcf.table.indels.clean.tmp out.vcf.table.indels.clean   
awk ' { for (i=6; i<=NF; i++) {
        if ($i ~ /|/) { 
          split($i, a, "|");
        if (a[1] == a[2]) $i=a[1];
           }
         }
       }; 
       {print $0} ' out.vcf.table.indels.clean > out.vcf.table.indels.clean.tmp		   
mv out.vcf.table.indels.clean.tmp out.vcf.table.indels.clean   

#remove any remaining mixed alleles
grep -v '|' out.vcf.table.indels.clean > out.vcf.table.indels.clean.tmp
mv out.vcf.table.indels.clean.tmp out.vcf.table.indels.clean  
grep -v '/' out.vcf.table.indels.clean > out.vcf.table.indels.clean.tmp
mv out.vcf.table.indels.clean.tmp out.vcf.table.indels.clean  
sed -i 's/ /\t/g' out.vcf.table.indels.clean

taxa=$(head -n1 out.vcf.table | cut -f3,6- | sed 's/.GT//g')
ntaxa=$(awk '{print NF-4; exit }' out.vcf.table)
nchar=$(cat out.vcf.table.indels.clean | wc -l)
awk '{print $1,$2}' out.vcf.table.indels.clean | sed 's/ /_/g' > indel.location
cut -f3,6- out.vcf.table.indels.clean > grid.indel.nucleotide
grid=$(paste indel.location grid.indel.nucleotide)
echo -e "\n#nexus\nbegin data;\ndimensions ntax=$ntaxa nchar=$nchar;\nformat symbols=\"AGCT\" gap=. transpose;\ntaxlabels $taxa;\nmatrix\n$grid\n;\nend;" > indel_matrix.nex


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

	
snpEff eff -no-downstream -no-intergenic -ud 100 -formatEff -v ${variant_genome_path} out.vcf > out.annotated.vcf
	
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

###########################################################################
# Creates the SNP and INDEL matrix for PAUP (or other phylogenetic programs)
###########################################################################

if [ ! -s "indel_matrix.nex" ]; then
    echo -e "\nScript must be supplied with indel_matrix.nex. Please check the directory and analysis and run again\n"
else
	echo -e "Found indel matrix\n"
fi

if [ ! -s "Ortho_SNP_matrix.nex" ]; then
	echo -e "Script must be supplied with Ortho_SNP_matrix.nex\n"
	echo -e "Please check the directory and analysis and run again\n"
else
	echo -e "Found Ortho_SNP_matrix.nex\n"
fi

#chomp files into just SNPs and positions
tail -n +8 Ortho_SNP_matrix.nex | head -n -2 > SNP_matrix_tmp
awk ' { for (i=3; i<=NF; i++) {if ($i == $2) $i=0; else $i=1}};  {print $0} ' SNP_matrix_tmp > SNP01.tmp

#indels
tail -n +8 indel_matrix.nex | head -n -2 > indel_matrix_tmp
awk ' { for (i=3; i<=NF; i++) {if ($i == $2) $i=0; else $i=1}};  {print $0} ' indel_matrix_tmp > indel01.tmp

cut -d " " -f 3- SNP01.tmp > SNP01.tmp2
awk '{ print $1 }' SNP01.tmp > SNP.loc
sed -i 's/$/ 0/g' SNP.loc 

cut -d " " -f 3- indel01.tmp > indel01.tmp2
awk '{ print $1 }' indel01.tmp > indel.loc
sed -i 's/$/ 0/g' indel.loc
paste -d ' ' indel.loc indel01.tmp2 > indel.mrg
paste -d ' ' SNP.loc SNP01.tmp2 > SNP.mrg

#cat files and create header
cat SNP.mrg indel.mrg > SNPindel.mrg
x=`cat SNPindel.mrg | wc -l`
y=`head -n 1 SNPindel.mrg | awk '{print NF}'`
z=$((y - 1))
taxa=`tail -n +6 Ortho_SNP_matrix.nex | head -n 1 |cut -d ' ' -f 2-`
grid=`cat SNPindel.mrg`
echo -e "\n#nexus\nbegin data;\ndimensions ntax=$z nchar=$x;\nformat symbols=\"01\" gap=. datatype=standard transpose;\ntaxlabels $taxa\nmatrix\n$grid\n;\nend;" > indel_SNP_matrix.nex

echo "SPANDx has finished"

exit 0
