#!/bin/bash

## This script will take two input files - inGroup outGroup - and construct the PED file for plink
#variable one should be inGroup
#variable number 2 should be outGroup
#
## The ped file has the following format
#column 1 = Family ID
#column 2 = Individual ID
#column 3 = Paternal ID
#column 4 = Maternal ID
#column 5 = sex (1 = male; 2 = female; other = unknown)
#column 6 = phenotype (0 = missing; 1 = unaffected; 2 = affected)
#column 7+ are the genotypes
#
##the map file has the following format
# column 1 = chromosome
# column 2 = rs or SNP identifier
# column 3 = genetic position
# column 4 = base pair position
#
# The inGroup file should contain all the strain names (matching with sequence read names) of isolates with the "affected" phenotype
# The outGroup file should contain all the strain names of isolates know to be "unaffected"
# All other isolates in the sequence directory will be classed as "missing"
#
#
# Strain_name \t Strain_name \t 0 \t 0 \t A A \t G G \t etc
#
#
#
#
# TO DO list
# include PLINK input from annotated outputs to include all indel variation
# Also include a gene by gene analysis to reduce noise low datasets
# Need better reference handling

usage()
{
echo -e  "\nUSAGE: GeneratePlink.sh -r <ReferenceGenome> -i <inGroup.txt> -o <OutGroup.txt> -c <P/A cutoff> -a <Annotation: yes/no> -v <Variant file in snpEff format if annotation is switched on>\n"
}
help()
{
echo -e "\nThanks for using SPANDx!!\n"
usage
echo -e "Please include a list of ingroup and outgroup strains\n"
echo -e "ingroup and outgroup lists must be correctly formatted with a single genome on each line followed by a line break\n"
echo -e "If this script fails to produce the expected files this is the most likely reason\n\nIncorrectly formatted files will easily break this script!\n"
echo -e "This script is designed to be run in the SPANDx analysis directory and relies on SPANDx output files to correctly function\n"
echo -e "If indel analysis has been run across your genomes this script will also generate a plink input file for these variants\n"

}

#Define path to SPANDx install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi


source "$SCRIPTPATH"/SPANDx.config 


OPTSTRING="hi:o:c:r:a:v:"

declare SWITCH
BedCutoff=$(echo "0.9" | bc)
ref=none
annotate=no
variant_genome=
reference=1

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		
		i) inGroup="$OPTARG"
		  echo "In group file = $inGroup"
		   ;;
        
		o) outGroup="$OPTARG"
		   echo "Out group file = $outGroup"
		   ;;
		   
		 c) BedCutoff=$(echo "$OPTARG" | bc)
			echo -e "\nCutoff for P/A analysis is $BedCutoff\n"	
			;;
		 r) ref="$OPTARG"
		    echo "Reference strain is $ref"
			;;
			
		a) annotate="$OPTARG"   
		   if [ "$annotate" == yes -o "$annotate" == no ]; then
			echo "A gene list file will be produced for input into PLINK"
		   else
		     echo -e "\nIf the optional -a parameter is used it must be set to yes or no\n"
			 echo -e "By default -a is set to no and variants will not be annotated\n\n"
		     usage
			 exit 1
           fi	
			;;
			
		 v) variant_genome="$OPTARG"
		   echo "SnpEff will use $variant_genome as the annotated reference"
		   echo "Please make sure this name matches the name in SnpEff database. Also verify chromosome names are correct. These can be found easily at ftp://ftp.ncbi.nih.gov/genomes"
           ;;
		   
		\?) help
		    exit 1
		    ;;
		
		h) help
		   exit 1
		   ;;
		   
		*) echo "script error: un-handled argument"
           help
		   exit 1
		   ;;
		   
		
	esac
done

if [ ! $PBS_O_WORKDIR ]
    then
        PBS_O_WORKDIR="$PWD"
fi

cd $PBS_O_WORKDIR


if [ "$ref" == "none" ]; then
	echo -e "You have not provided a reference. Do you want to continue?\n"
	echo -e "Type 'Y' to continue or press enter to exit\n"
	read ref_test
	#echo $ref_test
	if [ "$ref_test" == "Y" -o "$ref_test" == "y" -o "$ref_test" == "yes" ]; then
		echo -e "Continuing\n\n"
		reference=0
	else
		exit 1
	fi	
fi
## file checks and program checks
if [ "$reference" == 1 ]; then
	if [ ! -s "$ref.fasta" ]; then
        echo -e "Couldn't locate reference file. Please make sure that reference file is in the analysis directory,\n you have specified the reference name correctly, and that the .fasta extension is not included\n"
        exit 1
    else
	    echo -e "Found FASTA file\n"
	fi

ref_blank_lines=`grep -c '^$' $ref.fasta`

	if [ "$ref_blank_lines" != 0 ]; then
	    echo -e "ERROR: Reference FASTA file is formatted incorrectly and must contain 0 blank lines. Blank lines will cause BWA and GATK to fail."
        exit 1
    else
	    echo -e "FASTA file looks to contain no blank lines. Good.\n"
	fi
fi

java_test=`command -v "$JAVA"`

log_eval()
{
  cd $1
  echo -e "\nIn $1\n"
  echo "Running: $2"
  eval "$2"
  status=$?

  if [ ! $status == 0 ]; then
    echo "Previous command returned error: $status"
    exit 1
  fi
}

if [ -d tmp ]; then
	rm -f $PBS_O_WORKDIR/tmp/*
fi

if [ ! -d tmp ]; then
	mkdir tmp
fi

if [ "$(echo $BedCutoff '>' 1 | bc -l)" -eq 1 ]; then
  echo -e "Cutoff value for P/A analysis must be between 0 and 1\n"
  echo -e "Please change to a value between 0 and 1 or use the default value of 0.9\n"
  exit 1
fi

if [ "$annotate" == yes ]; then
    if [ ! -f "$SNPEFF" ]; then
	        echo "ERROR: Generate PLINK requires SnpEff to function. Please make sure the correct path is specified in SPANDx.config"
		    exit 1
    fi

	
	
fi




if [ "$annotate" == yes ]; then
    grep "$variant_genome" "$SNPEFF_CONFIG" &> /dev/null
    status=$?
    if [ ! $status == 0 ]; then
        echo "SPANDx couldn't find the reference genome in the SnpEff config file" 
		echo "The name of the annotated reference genome specified with the -v switch must match a reference genome in the SnpEff database"
        echo "Does the SnpEff.config file contain the reference specified with the -v switch?"
		echo "Is the SnpEff.config file in the location specified by SPANDx.config?"
	    echo "If both of these parameters are correct please refer to the SnpEff manual for further details on the setup of SnpEff"
		exit 1
    else
        echo -e "SPANDx found the reference file in the SnpEff.config file\n" 
    fi
	
	#test to see if the chromosome names in the SnpEff database match those of the reference file
	
	#CHR_NAME=`$JAVA -jar $SNPEFF dump "$variant_genome" | grep -A1 'Chromosomes names' | tail -n1 | awk '{print $2}'|sed "s/'//g"`
	#REF_CHR=`head -n1 "$ref".fasta | sed 's/>//'`  
	#if [ "$CHR_NAME" == "$REF_CHR" ]; then
	#    echo -e "Chromosome names in the SnpEff database match the reference chromosome names, good\n"
	#else
	 #   echo -e "Chromosome names in the SnpEff database DON'T match the reference chromosome names.\n"
		#echo -e "Please change the names of the reference file to match those in the SnpEff database.\n"
		#echo -e "If you are unsure what these are, run: $JAVA -jar $SNPEFF dump $variant_genome\n"
		#echo -e "The first chromosome name is $CHR_NAME.\n\n"
		#exit 1
	#fi	
	
	
	if [ ! -d "$SNPEFF_DATA/$variant_genome" ]; then
	    echo -e "Downloading reference genome to SnpEff database\n"
		echo -e "If the program hangs here please check that the proxy settings are correct and the cluster has internet access\n"
		echo -e "If required SNPEff databases can be manually downloaded and added to the SPANDx pipeline\n"
		echo -e "Running the following command:"
		echo "$JAVA $JAVA_PROXY -jar $SNPEFF download -v $variant_genome"
        echo -e "In the following directory $PBS_O_WORKDIR\n"		
		$JAVA ${JAVA_PROXY} -jar $SNPEFF download -v $variant_genome
	else 
        echo -e "Annotated reference database has already been downloaded for SnpEff\n"
    fi	
	#create bed file of genes
	echo -e "Creating BED file of genes for input into PLINK\n"
	if [ ! -s ${ref}.genes ]; then 
	    java -jar $SNPEFF genes2bed "$variant_genome" > "$ref".genes
	fi
fi


if [ ! -s "$inGroup" ]; then
    usage
	echo -e "\nScript must be supplied with a specified list of inGroup strains\n"
	echo -e "Script cannot find inGroup or file contains no data. Please check and rerun\n"
	
	exit 1
	else
	echo -e "The following file has been specified as the inGroup list \n\n$inGroup\n"
fi

if [ ! -s "$outGroup" ]; then
	echo -e "Script must be supplied with a list of outGroup strains\n"
	echo -e "Script cannot find outGroup or file contains no data. Please check and rerun\n"
	exit 1
	else
	echo -e "The following file has been specified as the outGroup list \n\n$outGroup\n"
fi

#clear arrays

unset $inGroupArrayTmp
unset $outGroupArrayTmp
unset $snps_Array
unset $ped_array
unset $inGroupArray
unset $outGroupArray


#read genome list for arrays
#assign in group and out group in what will become column 6

inGroupArrayTmp=(`cat $PBS_O_WORKDIR/$inGroup`)
outGroupArrayTmp=(`cat $PBS_O_WORKDIR/$outGroup`)
inGroupArray=("${inGroupArrayTmp[@]/%/ 0 0 0 2}")
outGroupArray=("${outGroupArrayTmp[@]/%/ 0 0 0 1}")
ped_array=("${outGroupArray[@]}" "${inGroupArray[@]}")

#Look for CRLF line breaks in incorrectly formatted reference files and attempts to fix formatting

grep -U $'\015' "$inGroup" &> /dev/null
status=$?
if [ $status == 0 ]; then
  echo -e "File looks to be windows formatted with CRLF breaks. Attempting to convert\n"
   tr -d '\015' < "$inGroup" > ${inGroup}.tmp
   mv ${inGroup}.tmp "$inGroup"
   rm ${inGroup}.tmp
fi

grep -U $'\015' "$outGroup" &> /dev/null
status=$?
if [ $status == 0 ]; then
  echo -e "File looks to be windows formatted with CRLF breaks. Attempting to convert\n"
   tr -d '\015' < "$outGroup" > ${outGroup}.tmp
   mv ${outGroup}.tmp "$outGroup"
   rm ${outGroup}.tmp
fi


#Genome inclusion test
#ingroup
n1=${#inGroupArrayTmp[@]}

if [ $n1 == 0 ]; then
        echo -e "Program couldn't find any genomes in the in group file"
        echo -e "Please check that genomes have been listed in this file and the file is in the correct format"
    	exit 1
fi

for (( i=0; i<n1; i++ )); do
	grep -P ^${inGroupArrayTmp[i]}'\t' Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex &> /dev/null
	status=$?
    if [ ! $status == 0 ]; then
        echo "I couldn't find all in group strains in the comparative SNPs files" 
		echo "Please check that all genome names are spelt correctly and are included in the analysis"
		echo "The first offending genome I encountered was ${inGroupArrayTmp[i]} but there may be others"
		exit 1
	fi
	grep -w ${inGroupArrayTmp[i]} ${outGroup} &> /dev/null
	status=$?
	if [ $status == 0 ]; then
		echo "There appears to be duplicate genome names in the input files i.e. One genome was found in both the ingroup and outgroup file. Please correct and rerun"
		echo "The first duplicate name encountered was ${inGroupArrayTmp[i]}"
		exit 1
	fi
done


n2=${#outGroupArrayTmp[@]}

if [ $n2 == 0 ]; then
        echo -e "Program couldn't find any genomes in the in group file"
        echo -e "Please check that genomes have been listed in this file and the file is in the correct format"
    	exit 1
fi

for (( i=0; i<n2; i++ )); do
	grep -P ^${outGroupArrayTmp[i]}'\t' Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex &> /dev/null
	status=$?
    if [ ! $status == 0 ]; then
        echo "I couldn't find all out group strains in the comparative SNPs files" 
		echo "Please check that all genome names are spelt correctly and are included in the analysis"
		echo "The first offending genome I encountered was ${outGroupArrayTmp[i]} but there may be others"
		exit 1
	fi	
done

#creation of gene by gene list

#Functions

GenerateAllIndels () {
#This function will take both the annotated SNP and annotated indel outputs and convert them into a gene by gene input into PLINK
cd $PBS_O_WORKDIR/tmp

printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg

cp $PBS_O_WORKDIR/Outputs/Comparative/All_indels_annotated.txt All_indels_annotated.txt

#prints column id for first non-genome column

FinalColumn=`head -n 1 All_indels_annotated.txt | awk '{ for(i=1; i<NF; i++) { if($i ~ /Binary_code/) {print i }}}'`

FinalColumn=`echo $(($FinalColumn - 1))`

cut  -f2-$FinalColumn < All_indels_annotated.txt | tail -n +2 > 012

sed -i 's/?/0/g' 012
sed -i 's/\./0/g' 012

awk ' BEGIN {FS=" "} { for (i=2; i<=NF; i++) {if ($i != $1 && $i != 0) $i=1; if ($i == $1 ) $i=2 }}; {print $0} ' 012 | sed 's/\t/ /g' | cut -d ' ' -f 2- > 012.tmp 

head -n 1 All_indels_annotated.txt | cut -f 3-$FinalColumn > taxa

awk 'BEGIN {FS=" "; OFS=""}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s %s",arr[j,i],(j==NR?"":OFS));} print "";}}' taxa > taxa.mrg

awk 'BEGIN {FS=" "; OFS=""}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s %s",arr[j,i],(j==NR?"":OFS));} print "";}}' 012.tmp > 012.mrg

paste taxa.mrg 012.mrg > 012_taxa

#fetch genomes contained in the in and out group files from the 012 outputs

awk '{ print $1 }' genomes | while read f; do grep -w ^$f' ' 012_taxa > $f.indels; done;
for f in *.indels; do sed -i 's/ //g' $f; done
for f in *.indels; do awk '{print $2 }' $f > $f.tmp; done;

for f in *.indels.tmp; do sed -i 's/./& &\t/g' $f; done;
indels_Array=(`ls *.indels.tmp`)
cat ${indels_Array[@]} > indels.mrg
paste genomes.mrg indels.mrg > All_indels.ped

mv All_indels.ped $PBS_O_WORKDIR/Outputs/Comparative/All_indels.ped

awk '{ print $1 } ' All_indels_annotated.txt | tail -n +2 > locs
awk 'BEGIN { FS = "_"; OFS = "_"} {$NF="";print $0 }' locs | rev | sed 's/_//' | rev > chr_names
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste chr_names locs ident basepairs > temp.map
awk '{print $1, $4, $2, $5 }' temp.map > All_indels.map
mv All_indels.map $PBS_O_WORKDIR/Outputs/Comparative

#cleanup

rm $PBS_O_WORKDIR/tmp/*

}


GenerateAllSnps () {
cd $PBS_O_WORKDIR/tmp

cp $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs_annotated.txt All_SNPs_annotated.txt

printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg
if [ ! -s All_SNPs_annotated.trans ]; then
	awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' All_SNPs_annotated.txt > All_SNPs_annotated.trans
fi
awk '{ print $1 }' genomes | while read f; do grep -P ^$f' ' All_SNPs_annotated.trans > $f.SNPs; done; 
for f in *.SNPs; do cut -d ' ' -f 2- < $f > $f.tmp; done; 

for f in *.SNPs.tmp; do sed -i 's/ //g' $f; done;

for f in *.SNPs.tmp; do sed -i 's/./& &\t/g' $f; done;

snps_Array=(`ls *.SNPs.tmp`)
cat ${snps_Array[@]} > SNPs.mrg

#need to test

sed -i 's/?/0/g' SNPs.mrg
sed -i 's/\./0/g' SNPs.mrg

paste genomes.mrg SNPs.mrg > All_SNPs.ped


mv All_SNPs.ped $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs.ped


# the map file has the following format
# column 1 = chromosome
# column 2 = rs or SNP identifier
# column 3 = genetic position
# column 4 = base pair position

awk '{ print $1 } ' All_SNPs_annotated.txt | tail -n +2 > locs
awk 'BEGIN { FS = "_"; OFS = "_"} {$NF="";print $0 }' locs | rev | sed 's/_//' | rev > chr_names
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste chr_names locs ident basepairs > temp.map
awk '{print $1, $4, $2, $5 }' temp.map > All_SNPs.map
mv All_SNPs.map $PBS_O_WORKDIR/Outputs/Comparative

#cleanup

rm $PBS_O_WORKDIR/tmp/*
}

GeneratePedFile () {


#all commands below are in tmp directory

cd $PBS_O_WORKDIR/tmp

printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg

cp $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex Ortho_SNP_matrix_RAxML.nex

# The perl expression of grep !!should!! only grab genome explicitly matching your provided genome names. If something breaks look here first!! 


awk '{ print $1 }' genomes | while read f; do grep -P ^$f'\t' Ortho_SNP_matrix_RAxML.nex > $f.SNPs; done; 
for f in *.SNPs; do awk '{print $2 }' $f > $f.tmp; done;
for f in *.SNPs.tmp; do sed -i 's/./& &\t/g' $f; done;

snps_Array=(`ls *.SNPs.tmp`)
cat ${snps_Array[@]} > SNPs.mrg
paste genomes.mrg SNPs.mrg > Ortho_SNPs.ped


mv Ortho_SNPs.ped $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.ped

#cleanup

rm $PBS_O_WORKDIR/tmp/*


}

GenerateMapFile () {

cd $PBS_O_WORKDIR/tmp

cp $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex Ortho_SNP_matrix.nex
awk '{ print $1 } ' Ortho_SNP_matrix.nex | tail -n +8 | head -n -2 > locs
awk 'BEGIN { FS = "_"; OFS = "_"} {$NF="";print $0 }' locs | rev | sed 's/_//' | rev > chr_names
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste chr_names locs ident basepairs > temp.map
awk '{print $1, $4, $2, $5 }' temp.map > Ortho_SNPs.map # column 4 should be the bp positions

mv Ortho_SNPs.map $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.map

#cleanup

rm $PBS_O_WORKDIR/tmp/*

}

GenerateIndelPed () {

cd $PBS_O_WORKDIR/tmp

printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg

cp $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex indel_matrix.nex

cat indel_matrix.nex | tail -n +8 | head -n -2 | cut -d ' ' -f 2- > 012

awk '  { for (i=2; i<=NF; i++) {if ($i != $1) $i=1; if ($i == $1 ) $i=2 }}; {print $0} ' 012 | cut -d ' ' -f 2- > 012.tmp 

sed -n -e '6p' indel_matrix.nex | sed  's/^taxlabels //' | sed 's/;$//' > taxa

awk 'BEGIN {FS=" "; OFS=""}{for (i=2;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=2;i<=big;i++){for(j=1;j<=NR;j++){printf("%s %s",arr[j,i],(j==NR?"":OFS));} print "";}}' taxa > taxa.mrg

awk 'BEGIN {FS=" "; OFS=""}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s %s",arr[j,i],(j==NR?"":OFS));} print "";}}' 012.tmp > 012.mrg

paste taxa.mrg 012.mrg > 012_taxa

#fetch genomes contained in the in and out group files from the 012 outputs

awk '{ print $1 }' genomes | while read f; do grep -w ^$f' ' 012_taxa > $f.indels; done;
for f in *.indels; do sed -i 's/ //g' $f; done
for f in *.indels; do awk '{print $2 }' $f > $f.tmp; done;

for f in *.indels.tmp; do sed -i 's/./& &\t/g' $f; done;
indels_Array=(`ls *.indels.tmp`)
cat ${indels_Array[@]} > indels.mrg
paste genomes.mrg indels.mrg > Ortho_indels.ped

mv Ortho_indels.ped $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.ped

rm $PBS_O_WORKDIR/tmp/* 
}

GenerateIndelMap () {

cp $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex indel_matrix.nex
awk '{ print $1 } ' indel_matrix.nex | tail -n +8 | head -n -2 > locs
awk 'BEGIN { FS = "_"; OFS = "_"} {$NF="";print $0 }' locs | rev | sed 's/_//' | rev > chr_names
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste chr_names locs ident basepairs > temp.map
awk '{print $1, $4, $2, $5 }' temp.map > Ortho_indels.map 

mv Ortho_indels.map $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.map

#cleanup

rm locs indel_matrix.nex temp.map ident basepairs basepairstmp

}

GeneratePAPed () {

cd $PBS_O_WORKDIR/tmp


cp $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt Bedcov_merge.txt

#print from col 4 to end of file

cut -d $'\t' -f4- Bedcov_merge.txt > bed_tmp1


#transpose
awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' bed_tmp1 > bed_tmp1.trans

#replace the bed file with the specified cutoffs to create 0 1 matrix
awk -v cutoff="$BedCutoff" '  { for (i=2; i<=NF; i++) {if ($i == 1.0) $i=1; if ($i >= cutoff) $i=1; if ($i < cutoff) $i=2 }}; {print $0} ' bed_tmp1.trans > bed_tmp2.trans


#fetch genomes contained in the in and out group files from the 012 outputs


printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg
awk '{ print $1 }' genomes | while read f; do grep -P ^$f' ' bed_tmp2.trans > $f.PA; done;

#The above line may break at certain times depending on genome names. There may be difficulties with grep thinking that _ are the start of complete names.

for f in *.PA ; do cut -d ' ' -f2- $f > $f.tmp; done;

for f in *.PA.tmp; do sed -i 's/./& &\t/g' $f; done;
PA_Array=(`ls *.PA.tmp`)
cat ${PA_Array[@]} > PA.mrg
paste genomes.mrg PA.mrg > PA_matrix.ped


#mapfile

awk '{ print $1 } ' Bedcov_merge.txt | tail -n +2 > chr_PA
awk '{ print $1, $2, $3 } ' Bedcov_merge.txt | tail -n +2 > locs_PA
sed -i 's/ /_/g' locs_PA
awk '{ print $2, $3 } '  Bedcov_merge.txt  | tail -n +2 > locs_col3_4_PA
sed -i 's/ /_/g' locs_col3_4_PA
paste -d ' ' chr_PA locs_PA locs_col3_4_PA locs_col3_4_PA > PA_matrix.map


mv PA_matrix.ped $PBS_O_WORKDIR/Outputs/Comparative/PA_matrix.ped
mv PA_matrix.map $PBS_O_WORKDIR/Outputs/Comparative/PA_matrix.map

rm $PBS_O_WORKDIR/tmp/*

}


if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.ped ]; then
    echo -e "Generating orthologous SNP ped file\n"
    GeneratePedFile
fi

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.map ]; then
    echo -e "Generating orthologous SNP map file\n"
    GenerateMapFile
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.ped ]; then
    echo -e "Generating orthologous indel ped file\n"
	echo -e "WARNING Generate PLINK script currently doesn't support the inclusion of a reference file in the GWAS output for indels\n"
    GenerateIndelPed
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.map ]; then
    echo -e "Generating orthologous indel map file\n"
	echo -e "WARNING Generate PLINK script currently doesn't support the inclusion of a reference file in the GWAS output for indels\n"
	GenerateIndelMap
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/PA_matrix.ped ]; then
    echo -e "Generating PA ped and map file\n"
	echo -e "WARNING Generate PLINK script currently doesn't support the inclusion of a reference file in the GWAS output for PA analysis\n"
    GeneratePAPed
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs_annotated.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs.ped ]; then
	echo -e "Generating All SNPs ped and map file\n"
	GenerateAllSnps
fi	

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/All_indels_annotated.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/All_indels.ped ]; then
	echo -e "Generating All indels ped and map file\n"
	GenerateAllIndels
fi	


echo -e "-----------------------------\n"
echo -e "Script has finished running\n"
echo -e "Output files are in Outputs/Comparative/\n"
echo -e "Happy GWAing...\n"
echo -e "Please refer to the excellent PLINK manual for how best to interrogate these files\n"
echo -e "-----------------------------\n"

exit 0
