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
#
#
#
#

usage()
{
echo -e  "USAGE: GeneratePlink.sh -i <inGroup.txt> -o <OutGroup.txt> -c <P/A cutoff>"
}
help()
{
echo -e "\nThanks for using SPANDx!!\n"
usage
echo -e "Please include a list of ingroup and outgroup strains\n"
echo -e "ingroup and outgroup lists must be correctly encoded (UTF-8) with a single genome on each line followed by a line break\n"
echo -e "If this script fails to produce the expected files this is the most likely reason\n"
echo -e "This script is designed to be run in the SPANDx analysis directory and relies on SPANDx output files to correctly function\n"
echo -e "If indel analysis has been run across your genomes this script will also generate a plink input file for these variants\n"
echo -e "Future developments include a P/A data matrix constructed from the BEDcov outputs for variable gene GWAS\n"

}

#Define path to SPANDx install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi



OPTSTRING="hi:o:c:"

declare SWITCH
BedCutoff=$(echo "0.9" | bc)


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
			
		\?) usage
		    exit 1
		    ;;
		
		h) usage
		   help
		   exit 1
		   ;;
		   
		*) echo "script error: unhandled argument"
           usage
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
	rm -rf tmp
fi

if [ ! -d tmp ]; then
	mkdir tmp
fi

if [ "$(echo $BedCutoff '>' 1 | bc -l)" -eq 1 ]; then
  echo -e "Cutoff value for P/A analysis must be between 0 and 1\n"
  echo -e "Please change to a value between 0 and 1 or use the default value of 0.9\n"
  exit 1
fi



#test for bedcov and if present run P/A matrix creation


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


inGroupArrayTmp=(`cat $PBS_O_WORKDIR/$inGroup`)
outGroupArrayTmp=(`cat $PBS_O_WORKDIR/$outGroup`)
inGroupArray=("${inGroupArrayTmp[@]/%/ 0 0 0 2}")
outGroupArray=("${outGroupArrayTmp[@]/%/ 0 0 0 1}")
ped_array=("${outGroupArray[@]}" "${inGroupArray[@]}")


#read genome list for arrays
#assign in group and out group in what will become column 6



#Genome inclusion test
#ingroup
n1=${#inGroupArrayTmp[@]}

if [ $n1 == 0 ]; then
        echo -e "Program couldn't find any genomes in the in group file"
        echo -e "Please check that genomes have been listed in this file and the file is in the correct format"
    	exit 1
fi

for (( i=0; i<n1; i++ )); do
	grep -w -P ${inGroupArrayTmp[i]}'\t' Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex &> /dev/null
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
	grep -w -P ${outGroupArrayTmp[i]}'\t' Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex &> /dev/null
	status=$?
    if [ ! $status == 0 ]; then
        echo "I couldn't find all out group strains in the comparative SNPs files" 
		echo "Please check that all genome names are spelt correctly and are included in the analysis"
		echo "The first offending genome I encountered was ${outGroupArrayTmp[i]} but there may be others"
		exit 1
	fi	
done

GeneratePedFile ()
{


#all commands below are in tmp directory

cd $PBS_O_WORKDIR/tmp



printf "%s\n" "${ped_array[@]}" > temp.col

awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes

sort -d genomes > genomes.mrg

cp $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex Ortho_SNP_matrix_RAxML.nex

awk '{ print $1 }' genomes | while read f; do grep -w -P $f'\t' Ortho_SNP_matrix_RAxML.nex > $f.SNPs; done; # the grep -w makes sure to grab only the exact matches of the grep pattern rather than an inclusive pattern 
for f in *.SNPs; do awk '{print $2 }' $f > $f.tmp; done;
for f in *.SNPs.tmp; do sed -i 's/./& &\t/g' $f; done;

snps_Array=(`ls *.SNPs.tmp`)
cat ${snps_Array[@]} > SNPs.mrg
paste genomes.mrg SNPs.mrg > Ortho_SNPs.ped


mv Ortho_SNPs.ped $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.ped

#cleanup

rm genomes temp.col genomes.mrg 
rm *.SNPs.tmp
rm *.SNPs
rm SNPs.mrg
rm Ortho_SNP_matrix_RAxML.nex


}
#Part two of the script is creation of the map file

GenerateMapFile ()
{

cd $PBS_O_WORKDIR/tmp

cp $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex Ortho_SNP_matrix.nex
awk '{ print $1 } ' Ortho_SNP_matrix.nex | tail -n +8 | head -n -2 > locs
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste locs ident basepairs > temp.map
awk '{print $2, $3, $1, $4 }' temp.map > Ortho_SNPs.map # column 4 should be the rs IDs

#there is an error with the above awk command for some datasets. Either column 4 or column 5 should be printed for the rs IDs but this appears to change? Need to find a better way of doing this. The underscore in the sed command for creation of the locs file is the issue with certain formatting of chromosome names

#the above error should be fixed now

mv Ortho_SNPs.map $PBS_O_WORKDIR/Outputs/Comparative

#cleanup

rm locs Ortho_SNP_matrix.nex temp.map ident basepairstmp basepairs

}


GenerateIndelPed ()
{

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

awk '{ print $1 }' genomes | while read f; do grep -w $f' ' 012_taxa > $f.indels; done;
for f in *.indels; do sed -i 's/ //g' $f; done
for f in *.indels; do awk '{print $2 }' $f > $f.tmp; done;

for f in *.indels.tmp; do sed -i 's/./& &\t/g' $f; done;
indels_Array=(`ls *.indels.tmp`)
cat ${indels_Array[@]} > indels.mrg
paste genomes.mrg indels.mrg > Ortho_indels.ped

mv Ortho_indels.ped $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.ped

rm $PBS_O_WORKDIR/tmp/* 
}

GenerateIndelMap ()
{

cp $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex indel_matrix.nex
awk '{ print $1 } ' indel_matrix.nex | tail -n +8 | head -n -2 > locs
sed 's/^/rs/g' locs > ident
sed 's/_/ /g' locs > basepairstmp
awk '{print $NF }' basepairstmp > basepairs
sed -i 's/^/0 /g' locs
paste locs ident basepairs > temp.map
awk '{print $2, $3, $1, $4 }' temp.map > Ortho_indels.map 

mv Ortho_indels.map $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.map

#cleanup

rm locs indel_matrix.nex temp.map ident basepairs basepairstmp

}


GeneratePAPed()
{

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
awk '{ print $1 }' genomes | while read f; do grep -w -P $f' ' bed_tmp2.trans > $f.PA; done;

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





#tests for ortho SNPs ped and map fiels and runs above commands if it doesn't exist

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.ped ]; then
    echo -e "Generating SNP ped file\n"
    GeneratePedFile
fi

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNPs.map ]; then
    echo -e "Generating SNP map file\n"
    GenerateMapFile
fi


#tests for ortho indels ped and map fiels and runs above commands if it doesn't exist
if [ -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.ped ]; then
    echo -e "Generating indel ped file\n"
    GenerateIndelPed
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_indels.map ]; then
    echo -e "Generating indel map file\n"
	GenerateIndelMap
fi

if [ -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/PA_matrix.ped ]; then
    echo -e "Generating PA ped and map file\n"
    GeneratePAPed
fi






echo -e "-----------------------------\n"
echo -e "Script has finished running\n"
echo -e "Output files are in Outputs/Comparative/\n"
echo -e "Happy GWAing...\n"
echo -e "-----------------------------\n"

exit 0

#general commands for plink
#-- how to run plink

#test the input files
plink --file Ortho_SNPs --allow-no-sex --noweb
plink --file Ortho_indels --allow-no-sex --noweb
plink --file PA_matrix --allow-no-sex --noweb

##make a binary PED file
plink --file Ortho_SNPs --make-bed --out Ortho1 --allow-no-sex --noweb
plink --file Ortho_indels --make-bed --out OrthoIndels --allow-no-sex --noweb
plink --file PA_matrix --make-bed --out OrthoPA --allow-no-sex --noweb


#range of significance values adjusted for multiple testing
plink --bfile Ortho1 --assoc --adjust --out as2 --allow-no-sex --noweb
plink --bfile OrthoIndels --assoc --adjust --out asindels --allow-no-sex --noweb
plink --bfile OrthoPA --assoc --adjust --out asPA --allow-no-sex --noweb
