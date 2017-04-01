#generate PLINK files for Roary output


#script should be run as follows
#GeneratePLINK_Roary.sh -i ingroupfile.txt -o outgroupfile.txt -r Roary_input 

help() {
cat << _EOF_

The script should be run as follows
GeneratePLINK_Roary.sh -i ingroupfile.txt -o outgroupfile.txt -r Roary_output (assumed to be gene_presence_absence.csv by default)

This script will convert the output of Roary to PLINK format for microbial Genome Wide Association

_EOF_

}





##PATH=/usr/local/R_v3.2.2/lib64/R/bin/:$PATH #specific for Cheetah HPC 



if [ ! $PBS_O_WORKDIR ]; then
        PBS_O_WORKDIR="$PWD"
fi

if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi

OPTSTRING="hi:o:r:"

declare SWITCH
BedCutoff=$(echo "0.9" | bc)
ref=none
annotate=no
variant_genome=""
reference=1
ROARY_FILE=gene_presence_absence.csv

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		
		i) inGroup="$OPTARG"
		  echo "In group file = $inGroup"
		   ;;
        
		o) outGroup="$OPTARG"
		   echo "Out group file = $outGroup"
		   ;;
		   
		r) ROARY_FILE="$OPTARG"
			echo -e "Roary output file = $ROARY_FILE\n"	
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




if [ ! -s "$inGroup" ]; then
    help
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

if [ -s "$ROARY_FILE" ]; then
    echo -e "Found pangenome file\n"
	else
    echo -e "Couldn't find pangenome file. Exiting\n"
	exit 1
fi	

unset inGroupArrayTmp
unset outGroupArrayTmp
unset snps_Array
unset ped_array
unset inGroupArray
unset outGroupArray


inGroupArrayTmp=(`cat $inGroup`)
outGroupArrayTmp=(`cat $outGroup`)
inGroupArray=("${inGroupArrayTmp[@]/%/ 0 0 0 2}")
outGroupArray=("${outGroupArrayTmp[@]/%/ 0 0 0 1}")
ped_array=("${outGroupArray[@]}" "${inGroupArray[@]}")

  
#Look for CRLF line breaks and attempts to fix formatting
grep -U $'\015' "$inGroup" &> /dev/null
status=$?
if [ $status == 0 ]; then
  echo -e "File looks to be windows formatted with CRLF breaks. Attempting to convert\n"
   tr -d '\015' < "$inGroup" > ${inGroup}.tmp
   mv ${inGroup}.tmp "$inGroup"
fi

grep -U $'\015' "$outGroup" &> /dev/null
status=$?
if [ $status == 0 ]; then
  echo -e "File looks to be windows formatted with CRLF breaks. Attempting to convert\n"
   tr -d '\015' < "$outGroup" > ${outGroup}.tmp
   mv ${outGroup}.tmp "$outGroup"
fi
  
  
#Genome inclusion test
#ingroup

n1=${#inGroupArrayTmp[@]}

if [ $n1 == 0 ]; then
        echo -e "Program couldn't find any genomes in the in group file"
        echo -e "Please check that genomes have been listed in this file and the file is in the correct format"
    	exit 1
fi



n2=${#outGroupArrayTmp[@]}

if [ $n2 == 0 ]; then
        echo -e "Program couldn't find any genomes in the in group file"
        echo -e "Please check that genomes have been listed in this file and the file is in the correct format"
        exit 1
fi

  echo -e "Generating PLINK input files from Roary output\n"

## To do: because the Roary output includes the core genome we should remove all core regions present in every single strain.
## Removal of the core region will prevent the Bonferroni correction being overly harsh with our significance testing
##create input from Roary output
  for (( i=0; i<n1; i++ )); do
	head -n1 "$ROARY_FILE" | egrep "\""${inGroupArrayTmp[i]}"\"" &> /dev/null
	status=$?
    if [ ! $status == 0 ]; then
        echo "I couldn't find all in group strains in the comparative SNPs files" 
		echo "Please check that all genome names are spelt correctly and are included in the analysis"
		echo "The first offending genome I encountered was ${inGroupArrayTmp[i]} but there may be others"
		exit 1
	fi
  done
  echo -e "Found all ingroup strains\n"  
  for (( i=0; i<n2; i++ )); do
      head -n1 "$ROARY_FILE" | grep "\""${outGroupArrayTmp[i]}"\"" &> /dev/null
      status=$?
      if [ ! $status == 0 ]; then
        echo "I couldn't find all out group strains in the comparative SNPs files" 
        echo "Please check that all genome names are spelt correctly and are included in the analysis"
        echo "The first offending genome I encountered was ${outGroupArrayTmp[i]} but there may be others"
        exit 1
	   fi     
  done
  echo -e "Found all outgroup strains\n"
  
  grep -w ${inGroupArrayTmp[i]} ${outGroup} &> /dev/null
  status=$?
  if [ $status == 0 ]; then
	echo "There appears to be duplicate genome names in the input files i.e. One genome was found in both the ingroup and outgroup file. Please correct and rerun"
	echo "The first duplicate name encountered was ${inGroupArrayTmp[i]}"
    exit 1
	else
	echo -e "No duplicates in the ingroup or outgroup files, good\n"
  fi
   
  echo -e "Creating R script\n"
 
  cat << '_EOF_' >> Recode_Roary.R
#!/usr/bin/env Rscript 

getArgs <- function(verbose=FALSE, defaults=NULL) {
  myargs <- gsub("^--","",commandArgs(TRUE))
  setopts <- !grepl("=",myargs)
  if(any(setopts))
    myargs[setopts] <- paste(myargs[setopts],"=notset",sep="")
  myargs.list <- strsplit(myargs,"=")
  myargs <- lapply(myargs.list,"[[",2 )
  names(myargs) <- lapply(myargs.list, "[[", 1)

  ## logicals
  if(any(setopts))
    myargs[setopts] <- TRUE

  ## defaults
  if(!is.null(defaults)) {
    defs.needed <- setdiff(names(defaults), names(myargs))
    if(length(defs.needed)) {
      myargs[ defs.needed ] <- defaults[ defs.needed ]
    }
  }

  ## verbage
  if(verbose) {
    cat("read",length(myargs),"named args:\n")
    print(myargs)
  }
  myargs
}

args <- getArgs()
## Default setting when no arguments passed
matrix <- args$matrix
output <- args$output

print(matrix)
print(output)



#load roary into R
library("dplyr")

my.data <- read.csv(matrix, stringsAsFactors=F, check.names = FALSE, na.strings="")

strain.sub <- select(my.data, -c(1:14))


#recode
strain.sub[is.na(strain.sub)] <- 2 
strain.sub[strain.sub != 2] <- 1

#transpose
strain.sub.t <- t(strain.sub)


#Would be good to drop all columns that have only "1" at this stage
#export

write.table(strain.sub.t, output, sep="\t", quote=FALSE)


_EOF_
  
  chmod +x $PBS_O_WORKDIR/Recode_Roary.R
  
  echo -e "Rscript created. Recoding Roary matrix\n"

  echo "$PBS_O_WORKDIR/Recode_Roary.R --no-save --no-restore --args matrix=$ROARY_FILE output=$PBS_O_WORKDIR/strain.data.trans.txt"
  
  $PBS_O_WORKDIR/Recode_Roary.R --no-save --no-restore --args matrix=$ROARY_FILE output=$PBS_O_WORKDIR/strain.data.trans.txt
  
  printf "%s\n" "${ped_array[@]}" > temp.col
  awk '{print $1, $1, $2, $3, $4, $5 }' temp.col > genomes
  sort -d genomes > genomes.mrg
  awk '{ print $1 }' genomes | while read f; do grep -P "^$f\t" strain.data.trans.txt > $f.PA; done;
  for f in *.PA ; do cut -d$'\t' -f2- $f > $f.tmp; done;
  for f in *.PA.tmp; do sed -i 's/./& &\t/g' $f; done;
  PA_Array=(`ls *.PA.tmp`)
  cat ${PA_Array[@]} > PA.mrg
  paste genomes.mrg PA.mrg > Roary_matrix.ped
  
  #mapfile
  awk 'BEGIN {FS=","} { print $1 }' "$ROARY_FILE"  | tail -n +2 | sed 's/ //g' | sed 's/"//g' > genes_PA
  lines=`wc -l genes_PA | awk '{ print $1 }'`
  #echo $lines

  for (( i=0; i<lines; i++ )); do
     pos=`echo "$i * 100" | bc`
     echo -e "1\t$pos" >> chromosome
  done

  paste -d ' ' genes_PA chromosome > tmp

  awk 'BEGIN {FS=OFS=" "} {print $2, $1, $3, $3}' tmp > Roary_matrix.map
  
  rm *.tmp *.PA tmp chromosome genes_PA PA.mrg temp.col genomes.mrg genomes strain.data.trans.txt Recode_Roary.R
  
 # $plink --file Roary_matrix --make-bed --out Roary --allow-no-sex --allow-extra-chr
  #$plink --bfile Roary --assoc --adjust --out asPresenceAbsence --allow-no-sex --allow-extra-chr
  
  
  #skipped pdf file creation for Roary pangenome. Need to correctly handle location of region in the creation of the BP file
  
exit 0
