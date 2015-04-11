#!/bin/bash
#$ -S /bin/bash
#PBS -S /bin/sh
#$ -cwd

#########################################################################
# The following script will combine all vcf files in your GATK directory into a master VCF file
# Combined snp calls from the master.vcf are interrogated across all bam files at these locations to verify snp identity
# Finally the clean vcf files are concatenated and converted into a matrix for phylogeny programs
# Version 1.4
# Written by D. Sarovich
# derek.sarovich@menzies.edu.au
#
# Version history
# 1.2-1.3 - added SGE bash interpreter line to the header
# Version 1.1-1.2 
# Added post processing of SNP matrix to remove lines containing 1 and -1 and combined individual samples and IDs into a single matrix
# Some additional file cleanup 
## V1.4
# added SGE header
#########################################################################



#source variables
source "$SCRIPTPATH"/SPANDx.config

if [ ! $PBS_O_WORKDIR ]; then
        PBS_O_WORKDIR="$seq_path"
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
#########################################################################
## checks and creates directory structure
#########################################################################

if [ ! -d "$seq_path/Phylo" ]; then
    mkdir $seq_path/Phylo
  else
    echo -e "Phylo directory already exists \n"
fi
if [ ! -d "$seq_path/Phylo/snps" ]; then
    mkdir $seq_path/Phylo/snps
  else
    echo -e "Phylo/snps directory already exists\n"
fi  
if [ ! -d "$seq_path/Phylo/bams" ]; then
    mkdir $seq_path/Phylo/bams
  else
    echo -e "Phylo/bams directory already exists\n"
fi  
if [ ! -d "$seq_path/Phylo/out" ]; then
    mkdir $seq_path/Phylo/out
  else
    echo -e "Phylo/out directory already exists\n"
fi  

if [ ! -d "$seq_path/Phylo/indels" ]; then
    mkdir $seq_path/Phylo/indels
  else
    echo -e "Phylo/indels/out directory already exists\n"
fi

#########################################################################
#links and rename snps PASS files into snps directory
#########################################################################

if [ "$(ls -A $PBS_O_WORKDIR/Phylo/snps)" ]; then
   echo "snps directory not empty"
   echo "skipping linking and renaming of vcf files"
   echo -e "if the snp vcf files are not in snps directory, make sure the directory is empty and restart\n"
  else
   log_eval $PBS_O_WORKDIR/Phylo/snps "ls ../../Outputs/SNPs_indels_PASS/*.snps.PASS.vcf|while read f; do ln -s \$f; done;"
   log_eval $PBS_O_WORKDIR/Phylo/snps "for f in *.snps.PASS.vcf; do mv \$f \${f//.snps.PASS.vcf/.vcf}; done;"
fi
if [ "$(ls -A $PBS_O_WORKDIR/Phylo/indels)" ]; then
   echo "indels directory not empty"
   echo "skipping linking and renaming of indel vcf files"
   echo -e "if the indel vcf files are not in indels directory, make sure the directory is empty and restart\n"
  else
   log_eval $PBS_O_WORKDIR/Phylo/indels "ls ../../Outputs/SNPs_indels_PASS/*.indels.PASS.vcf|while read f; do ln -s \$f; done;"
   log_eval $PBS_O_WORKDIR/Phylo/indels "for f in *.indels.PASS.vcf; do mv \$f \${f//.indels.PASS.vcf/.vcf}; done;"
fi
if [ "$(ls -A $PBS_O_WORKDIR/Phylo/bams)" ]; then
   echo "bams directory not empty"
   echo "skipping linking and renaming of bam and bai files"
   echo -e "if the bam and bai files are not in bams directory, make sure the directory is empty and restart\n"
  else
   log_eval $PBS_O_WORKDIR/Phylo/bams "ls ../../*/unique/*.realigned* |while read f; do ln -s \$f; done;"
   log_eval $PBS_O_WORKDIR/Phylo/bams "for f in *.bam; do mv \$f \${f//.realigned.bam/.bam}; done;"
   log_eval $PBS_O_WORKDIR/Phylo/bams "for f in *.bai; do mv \$f \${f//.realigned.bai/.bai}; done;"
fi
if [ ! -d "$seq_path/Phylo/indels/out" ]; then
    mkdir $seq_path/Phylo/indels/out
  else
    echo -e "Phylo/indels directory already exists\n"
fi 
#########################################################################
#Checks for the master vcf file and creates it
#If merge indels is switched on a master.vcf will be created for the indels as well
#########################################################################

if [ ! -s $PBS_O_WORKDIR/Phylo/out/master.vcf ]; then
    array=($(find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "))
    array2=("${array[@]/.vcf/,VCF}")
    array3=("${array2[@]/#/-V:}")
    n=${#array3[@]}
    for (( i=0; i<n; i++ )); do input[i]=${array3[i]}" "${array[i]}; done;
    log_eval $PBS_O_WORKDIR/Phylo/snps "$JAVA $SET_VAR $GATK -T CombineVariants -R $PBS_O_WORKDIR/${ref}.fasta ${input[*]} --sites_only -o $PBS_O_WORKDIR/Phylo/out/master.vcf"
  else
    echo -e "The merged snp calls have already been combined into the master.vcf file\n"
fi

if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/master.vcf -a "$indel_merge" == yes ]; then
    array=($(find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "))
    array2=("${array[@]/.vcf/,VCF}")
    array3=("${array2[@]/#/-V:}")
    n=${#array3[@]}
    for (( i=0; i<n; i++ )); do input[i]=${array3[i]}" "${array[i]}; done;
    log_eval $PBS_O_WORKDIR/Phylo/indels "$JAVA $SET_VAR $GATK -T CombineVariants -R $PBS_O_WORKDIR/${ref}.fasta ${input[*]} --sites_only -o $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf"
fi

sleep 20
 
exit 0
