#!/bin/bash

#####################################################
# 
# Thanks for using SPANDx!!
#
# USAGE: SPANDx.sh <parameters, required> -r <reference, without .fasta extension> [parameters, optional] -o [organism] -m [generate SNP matrix yes/no] -i [generate indel matrix yes/no] 
# -a [include annotation yes/no] -v [Variant genome file. Name must match the SnpEff database] -s [Specify read prefix to run single strain or none to construct a SNP matrix from a previous analysis ] -t [Sequencing technology used Illumina/Illumina_old/454/PGM] 
# -p [Pairing of reads PE/SE] -w [Window size in base pairs for BEDcoverage module]
#
# SPANDx by default expects reads to be paired end, in the following format: STRAIN_1_sequence.fastq.gz for the first pair and STRAIN_2_sequence.fastq.gz for the second pair. 
# Reads not in this format will be ignored.
# If your data are not paired, you must set the -p parameter to SE to denote unpaired reads. By default -p is set to PE.
#
# SPANDx expects at least a reference file in FASTA format. 
# For compatibility with all programs in SPANDx, FASTA files should conform to the specifications listed here: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml
# Note that the use of nucleotides other than A, C, G, or T is not supported by certain programs in SPANDx and these should not be used in reference FASTA files. 
# In addition, Picard, GATK and SAMtools handle spaces within contig names differently. Therefore, please avoid using spaces or special characters (e.g. $/*) in contig names.
# 
# By default all read files present in the current working directory will be processed. 
# Sequence files within current directory will be aligned against the reference using BWA, SNPs and indels will be called with GATK and a SNP 
# matrix will be generated with GATK and VCFTools
# 
# Optionally to run the SNP matrix without processing any sequence files set -s to none and -m to yes. If -s is set to none SPANDx will skip all sequence files in the current directory and will not perform the alignment or variant calls. Instead SPANDx will merge all VCF files contained in $PWD/Phylo/snps, interrogate those calls in the alignment files contained within $PWD/Phylo/bams and output a SNP matrix in $PWD/Phylo/out. Before running this module please check that all VCF files located within $PWD/Phylo/snps match the alignment files within $PWD/Phylo/bams
# e.g. SPANDx.sh -r K96243 -m yes -s none
#
# Written by Derek Sarovich and Erin Price - Menzies School of Health Research, Darwin, Australia
# Please send bug reports to mshr.bioinformatics@gmail.com
# If you find SPANDx useful and use it in published work please cite - SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618"
# Version 2.2
# 2.0-2.1 Added SGE job handling
# 2.1-2.2 Added SLURM job handling
#
#################################################################
usage()
{
echo -e  "USAGE: SPANDx.sh <parameters, required> -r <reference, without .fasta extension> [parameters, optional] -o [organism] -m [generate SNP matrix yes/no] -i [generate indel matrix yes/no] -a [include annotation yes/no] -v [Variant genome file. Name must match the SnpEff database] -s [Specify read prefix to run single strain or none to just construct SNP matrix] -t [Sequencing technology used Illumina/Illumina_old/454/PGM] -p [Pairing of reads PE/SE] -w [BEDTools window size in base pairs]"
}
help()
{
echo -e "\nThanks for using SPANDx!!\n"
usage
echo -e "The only essential option is to specify a reference. All other parameters are optional. By default the program will process all Next-Generation"
echo -e "sequencing reads within the present working directory.\n"
echo -e "SPANDx by default expects reads to be paired end, in the following format: STRAIN_1_sequence.fastq.gz for the first pair and STRAIN_2_sequence.fastq.gz for the second pair."
echo -e "Reads not in this format will be ignored"
echo -e "If your data are not paired, you must set the -p parameter to yes to denote unpaired reads. By default -p is set to PE.\n"
echo -e "SPANDx expects a reference file in FASTA format."
echo -e "Although SPANDx will accept FASTA files with very loose formatting specifications, for compatibility"
echo -e "with all programs FASTA files should conform to the specifications listed here: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml"
echo -e "NB. Some programs used within SPANDx do not allow the use of IUPAC codes and these should not be used in reference FASTA files\n"
echo -e "By default all read files present in the current working"
echo -e "directory will be processed. Sequences will be aligned against the reference using BWA, SNPs and indels will be called with GATK and a SNP"
echo -e "matrix will be generated with GATK and VCFTools\n"
echo -e "e.g. SPANDx.sh -r K96243 -o Bpseudo -m yes"
echo -e "Please send bug reports to mshr.bioinformatics@gmail.com\n\n"
echo -e "If you use SPANDx in published work please cite - SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618"
}

#Define path to SPANDx install
SCRIPTPATH="$( cd "$(dirname "$0")" ; pwd -P )"

if  [ $# -lt 1 ]
    then
	    usage
		exit 1
fi

# source dependencies
source "$SCRIPTPATH"/SPANDx.config 
source "$SCRIPTPATH"/scheduler.config

#declare variables
declare -rx SCRIPT=${0##*/}


OPTSTRING="hr:o:d:m:a:v:s:t:p:w:i:"

declare SWITCH
declare ref
org=haploid
seq_directory="$PWD"
matrix=no
annotate=no
strain=all
tech=Illumina
pairing=PE
variant_genome=
window=1000
indel_merge=no

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do 
		case $SWITCH in
		
		r) ref="$OPTARG"
		   echo "Reference = $ref"
		   ;;
        
		o) org="$OPTARG"
		   echo "Organism = $org"
		   ;;
		
		d) seq_directory="$OPTARG"
		   echo "SPANDx expects read files to be in $seq_directory"
		   ;;
		   
        m) matrix="$OPTARG"
           if [ "$matrix" == yes -o "$matrix" == no ]; then 
           echo "SNP matrix will be generated? $matrix"
		   else
		       echo -e "\nIf the optional -m parameter is used it must be set to yes or no\n"
			   echo -e "By default -m is set to no and SPANDx will not generate a SNP matrix\n"
		       usage
			   exit 1
           fi		
           ;;  
		   
		a) annotate="$OPTARG"   
		   if [ "$annotate" == yes -o "$annotate" == no ]; then
		   echo "VCF files will be annotated with SnpEff? $annotate"
		   else
		       echo -e "\nIf the optional -a parameter is used it must be set to yes or no\n"
			   echo -e "By default -a is set to no and variants will not be annotated\n\n"
		       usage
			   exit 1
           fi		
           ;; 
		
		s) strain="$OPTARG"
		   echo "Only strain $strain will be processed. By default all strains within current directory will be processed"
		   ;;
		
		t) tech="$OPTARG"
			if [ "$tech" == Illumina -o "$tech" == Illumina_old -o "$tech" == 454 -o "$tech" == PGM ]; then
			    echo -e "Technology used is $tech\n"
			else
			   echo -e "If the optional parameter -t is used it must be set to Illumina, Illumina_old, 454 or PGM.\n"
			   echo -e "By default -t is set to Illumina and SPANDx will assume your sequence data is in Illumina format with quality encoding of 1.8+\n"
			   usage
			   exit 1
			fi
			;;
			
		p) pairing="$OPTARG"
		   if [ "$pairing" == PE -o "$pairing" == SE ]; then 
               echo "The pairing of your reads has been indicated as $pairing"
		   else
		       echo -e "\nIf the optional -p parameter is used it must be set to PE or SE\n"
			   echo -e "By default -p is set to PE and SPANDx assumes your reads are paired end\n"
		       usage
			   exit 1
		   fi
		   ;;
		   
		w) window="$OPTARG"
		   echo "BEDTools BEDCoverage window size has been set to $window base pairs. Default is 1000bp"
		   ;;
		
		v) variant_genome="$OPTARG"
		   echo "SnpEff will use $variant_genome as the annotated reference"
		   echo "Please make sure this name matches the name in SnpEff database. Also verify chromosome names are correct. These can be found easily at ftp://ftp.ncbi.nih.gov/genomes"
           ;;
		
        i) indel_merge="$OPTARG"
		   if [ "$indel_merge" == yes -o "$indel_merge" == no ]; then
               echo -e "Indels will be merged and checked across genomes = $indel_merge\n"  
		   else
		       echo -e "Indel merge must be set to yes or no. Please refer to the manual for more details\n"
			   exit 1
			   fi
           ;;
		   
		\?) usage
		    exit 1
		    ;;
		
		h) usage
		   exit 1
		   ;;
		   
		*) echo "script error: unhandled argument"
           exit 1
		   usage
		   ;;
		   
		
	esac
done
  
echo -e "\nThe following parameters will be used\n"
echo -e "-------------------------------------\n"
echo "Organism = $org"
echo "Output directory and directory containing sequence files = $seq_directory"
echo "SNP matrix will be created? = $matrix"
echo "Genomes will be annotated? = $annotate"
echo "Strain(s) to be processed = $strain"
echo "Sequence technology used = $tech"
echo "Pairing of reads = $pairing"
echo "Variant genome specified for SnpEff = $variant_genome"
echo "Window size for BedTools = $window"
echo "Indels will be merged and corrected = $indel_merge"
echo -e "-------------------------------------\n\n"

ref_index=${seq_directory}/${ref}.bwt #index file created with BWA
REF_INDEX_FILE=${seq_directory}/${ref}.fasta.fai #index created with SAMtools
REF_BED=${seq_directory}/${ref}.bed
REF_DICT=${seq_directory}/${ref}.dict #Dictionary file created with Picard tools

if [ ! $PBS_O_WORKDIR ]; then
        PBS_O_WORKDIR="$seq_directory"
fi

cd $PBS_O_WORKDIR

## file checks and program checks

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

## Test for dependencies required by SPANDx


## Test for dependencies required by SPANDx
bwa_test=`command -v "$BWA"`
samtools_test=`command -v "$SAMTOOLS"`
bedtools_test=`command -v "$BEDTOOLS"`
vcftools_test=`command -v "$VCFTOOLS"`
vcfmerge_test=`command -v "$VCFMERGE"`
bgzip_test=`command -v "$BGZIP"`
tabix_test=`command -v "$TABIX"`
java_test=`command -v "$JAVA"`

if [ -z "$bwa_test" ]; then
	    echo "ERROR: SPANDx requires BWA to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$samtools_test" ]; then
	    echo "ERROR: SPANDx requires SAMtools to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ ! -f "$GATK" ]; then
	    echo "ERROR: SPANDx requires the Genome Analysis toolkit and java to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$bedtools_test" ]; then
	    echo "ERROR: SPANDx requires  BEDtools to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ ! -f "$PICARD" ]; then
	    echo "ERROR: SPANDx requires Picard to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
#if [ ! -f "$ADDORREPLACEREADGROUPS" ]; then
#	    echo "ERROR: SPANDx requires Picard to function. Please make sure the correct path is specified in SPANDx.config"
#		exit 1
#fi
#if [ ! -f "$BUILDBAMINDEX" ]; then
#	    echo "ERROR: SPANDx requires Picard to function. Please make sure the correct path is specified in SPANDx.config"
#		exit 1
#fi
#if [ ! -f "$CREATEDICT" ]; then
#	    echo "ERROR: SPANDx requires Picard to function. Please make sure the correct path is specified in SPANDx.config"
#		exit 1
#fi
if [ "$annotate" == yes ]; then
    if [ ! -f "$SNPEFF" ]; then
	        echo "ERROR: SPANDx requires SnpEff to function. Please make sure the correct path is specified in SPANDx.config"
		    exit 1
    fi
fi
if [ -z "$vcftools_test" ]; then
	    echo "ERROR: SPANDx requires VCFtools to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$vcfmerge_test" ]; then
	    echo "ERROR: SPANDx requires vcf-merge to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$bgzip_test" ]; then
	    echo "ERROR: SPANDx requires bgzip to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$tabix_test" ]; then
	    echo "ERROR: SPANDx requires tabix to function. Please make sure the correct path is specified in SPANDx.config"
		exit 1
fi
if [ -z "$java_test" ]; then
	    echo "ERROR: SPANDx requires java. Please make sure java is available on your system"
		exit 1
fi

# Sequence Technology check for pairing

if [ "$tech" == 454 -o "$tech" == PGM ]; then
    if [ "$pairing" == PE ]; then
	    echo -e "ERROR: SPANDx does not currently support paired end PGM or 454 data. If your PGM or 454 data is single end please change the -p switch to SE\n"
		exit 1
    fi
fi

### Handling and checks for read files
if [ "$strain" == all ]; then
    sequences_tmp=(`find $PBS_O_WORKDIR/*_1_sequence.fastq.gz -printf "%f "`)
    sequences=("${sequences_tmp[@]/_1_sequence.fastq.gz/}")
    n=${#sequences[@]}
    if [ $n == 0 ]; then
        echo -e "Program couldn't find any sequence files to process"
        echo -e "Sequences must be in the format STRAIN_1_sequence.fastq.gz STRAIN_2_sequence.fastq.gz for paired end reads"
    	echo -e "and STRAIN_1_sequence.fastq.gz for single end data\n"
	    exit 1
    else
        echo -e "Sequences have been loaded into SPANDx\n"
    fi
fi

## check for read pairing and correct notation #need to test
if [ "$pairing" == PE -a "$strain" == all ]; then
	sequences_tmp2=(`find $PBS_O_WORKDIR/*_2_sequence.fastq.gz -printf "%f "`)
    sequences2=("${sequences_tmp2[@]/_2_sequence.fastq.gz/}")
    n2=${#sequences2[@]}
    if [ $n != $n2 ]; then
	    echo "Number of forward reads don't match number of reverse reads. Please check that for running in PE mode all read files have correctly named pairs"
		exit 1
	fi
	for (( i=0; i<n; i++ )); do
	    if [ ${sequences[$i]} != ${sequences2[$i]} ]; then
            echo "Names of forward reads don't match names of reverse reads. Please check that for running in PE mode all read files have correctly named pairs"
			echo -e "The offending read pair is ${sequences[$i]} and ${sequences2[$i]}\n"
			#echo "Forward read names are ${sequences[@]}"
			#echo "Reverse read names are ${sequences2[@]}"
            exit 1
        fi
    done;
fi

## error checking for strain variable. Makes sure that if a single strain is specified that the read data is present

#to do 
#include if n == 1 and matrix == yes then exit
if [ "$strain" != all -a "$strain" != none ]; then
    sequences="$strain"
	echo -e "SPANDx will process the single strain $strain.\n"
	if [ "$matrix" == yes ]; then
	    matrix=no
		echo -e "SPANDx will not run SNP matrix with only one strain specified. If a SNP matrix is desired please run SPANDx with more strains or if strains have already been run link the SNP VCF and BAM files to the appropriate directories and set -s to none.\n"
	fi
	if [ "$pairing" == PE ]; then
	    if [ ! -s ${strain}_1_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_1_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
		if [ ! -s ${strain}_2_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_2_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
	fi
	if [ "$pairing" == SE ]; then
		if [ ! -s ${strain}_1_sequence.fastq.gz ]; then
		    echo -e "ERROR: SPANDx cannot find sequence file ${strain}_1_sequence.fastq.gz\n"
			echo -e "Please make sure the correct sequence reads are in the sequence directory\n"
			exit 1
        fi
	fi
fi	
	
#to do. Needs to have an n > 1 check before the test	
	
if [ "$matrix" == no -a "$indel_merge" == yes ]; then
    echo -e "Indel merge has been requested. As this cannot be determine without creation of a SNP matrix the matrix variable has been set to yes\n"
	matrix=yes
fi
	
	## create directory structure
	

if [ ! -d "BEDcov" ]; then
    mkdir $seq_directory/BEDcov 
fi
if [ ! -d "tmp" ]; then
	mkdir $seq_directory/tmp
fi
if [ ! -d "Outputs" ]; then
	mkdir $seq_directory/Outputs
fi
if [ ! -d "Outputs/SNPs_indels_PASS" ]; then
	mkdir $seq_directory/Outputs/SNPs_indels_PASS
fi
if [ ! -d "Outputs/SNPs_indels_FAIL" ]; then
	mkdir $seq_directory/Outputs/SNPs_indels_FAIL
fi
if [ ! -d "Outputs/Comparative" ]; then
	mkdir $seq_directory/Outputs/Comparative
fi

## checking variables for the annotation module

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
	if [ ! -d "$SNPEFF_DATA/$variant_genome" ]; then
	    echo -e "Downloading reference genome to SnpEff database\n"
		echo -e "If the program hangs here please check that the proxy settings are correct and the cluster has internet access\n"
		echo -e "If required SNPEff databases can be manually downloaded and added to the SPANDx pipeline\n"
		echo -e "Running the following command:"
        echo -e "In the following directory $PBS_O_WORKDIR\n"		
		echo "$JAVA $JAVA_PROXY -jar $SNPEFF download -v $variant_genome"
	    $JAVA ${JAVA_PROXY} -jar $SNPEFF download -v $variant_genome
	else 
        echo -e "Annotated reference database has already been downloaded for SnpEff\n"
    fi	
fi



# The following section will use qsub to run all of the SPANDx jobs specified in the command line
if [ "$SCHEDULER" == PBS ]; then

#clean batch system log files

if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi
if [ -s qsub_ids2.txt ]; then
    rm qsub_ids2.txt
fi
if [ -s qsub_array_ids.txt ]; then
    rm qsub_array_ids.txt
fi
if [ -s mastervcf_id.txt ]; then
    rm mastervcf_id.txt
fi
if [ -s clean_vcf_id.txt ]; then
	rm clean_vcf_id.txt
fi
if [ -s matrix_id.txt ]; then
    rm matrix_id.txt
fi


## Indexing of the reference with SAMTools and BWA
## creation of the BED file for the bedcoverage module
## creation of the GATK and picard reference dictionaries

  if [ ! -s "$ref_index" ]; then
	echo -e "Submitting qsub job for BWA reference indexing\n"
    cmd="$BWA index -a is -p ${seq_directory}/${ref} ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N index -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "index\t$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_INDEX_FILE" ]; then
    echo -e "Submitting qsub job for SAMtools reference indexing\n"
    cmd="$SAMTOOLS faidx ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N SAM_index -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "SAM_index\t$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_DICT" ]; then
    echo -e "Submitting qsub job for ${ref}.dict creation\n"
    cmd="$JAVA $SET_VAR $CREATEDICT R=${seq_directory}/${ref}.fasta O=$REF_DICT"
	qsub_id=`qsub -N PICARD_dict -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
	echo -e "PICARD_dict\t$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_BED" -a qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
	qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
	depend="-W depend=afterok${qsub_cat_ids}"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
   echo -e "BED_window\t$qsub_id" >> qsub_ids.txt	
  fi
  if [ ! -s "$REF_BED" -a ! qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "BED_window\t$qsub_id" >> qsub_ids.txt	
  fi

variants_single ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting qsub job for sequence alignment and variant calling for $sequences\n"
        var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N aln_sequences -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
        echo -e "aln_sequences\t$qsub_array_id" >> qsub_array_ids.txt
	fi
        
fi
if [ ! -s qsub_ids.txt ]; then
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N aln_$sequences -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
		echo -e "aln_$sequences\t$qsub_array_id" >> qsub_array_ids.txt
	fi
fi

}

## Variants is the main alignment and variant calling script contained with SPANDx

variants ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
                var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N aln_${sequences[$i]} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
                echo -e "aln_${sequences[$i]}\t$qsub_array_id" >> qsub_array_ids.txt
	        fi
        done
fi
if [ ! -s qsub_ids.txt ]; then
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    	    var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N aln_${sequences[$i]} -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
				echo -e "aln_${sequences[$i]}\t$qsub_array_id" >> qsub_array_ids.txt
	        fi
        done
fi

}


## Matrix is the main comparative genomics section of SPANDx that relies on the outputs from the variant function above
matrix ()			
{	
if [ -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> mastervcf_id.txt	
fi
if [ ! -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
	depend="-W depend=afterok${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
		qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
		echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done 
fi
## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes

if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    qsub_cat_ids=`cat clean_vcf_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
    depend="-W depend=afterok${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> matrix_id.txt
fi
}

## This function will generate a SNP matrix from the Phylo directory assuming all SNP and BAM files have already been linked into these directories
## qsub for indel_merge is currently untested

matrix_final ()			
{
if [ -s qsub_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'`
    depend="-W depend=afterok${qsub_cat_ids}"
	echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T $depend -v "$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> mastervcf_id.txt
fi
if [ ! -s qsub_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
	depend="-W depend=afterok${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes


if [ -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
	depend="-W depend=afterok${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "clean_vcf\t$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

#####################


if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    qsub_cat_ids=`cat clean_vcf_id.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
    depend="-W depend=afterok${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "Matrix_vcf\t$qsub_matrix_id" >> matrix_id.txt
fi
}

## This function takes the output of each bedcoverage assessment of the sequence alignment and merges them in a comparative matrix
## This function is run when $strain=all but not when a single strain is analysed i.e. $strain doesn't equal all
merge_BED ()
{
if [ -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f2 | sed -e 's/^/:/' | tr -d '\n'` 
    depend="-W depend=afterok${qsub_cat_ids}"
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T "$depend" -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "BEDcov_merge\t$qsub_BED_id" >> qsub_BED_id.txt
fi
if [ ! -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="ref=$ref,seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUTPUT -m $MAIL -M $ADDRESS -l ncpus=1,walltime=$WALL_T -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "BEDcov_merge\t$qsub_BED_id" >> qsub_BED_id.txt
fi
}

### These variable tests determine which of the above functions need to be run for the SPANDx pipeline

  if [ "$strain" == all ]; then
    variants
	merge_BED
  fi
  if [ "$strain" != all -a "$strain" != none ]; then
    variants_single
  fi
  if [ "$matrix" == yes -a "$strain" != none ]; then
    matrix
  fi
  if [ "$matrix" == yes -a "$strain" == none ]; then
    matrix_final
  fi
fi

if [ "$SCHEDULER" == SGE ]; then

#clean batch system log files

if [ -s qsub_ids.txt ]; then
    rm qsub_ids.txt
fi
if [ -s qsub_ids2.txt ]; then
    rm qsub_ids2.txt
fi
if [ -s qsub_array_ids.txt ]; then
    rm qsub_array_ids.txt
fi
if [ -s mastervcf_id.txt ]; then
    rm mastervcf_id.txt
fi
if [ -s clean_vcf_id.txt ]; then
	rm clean_vcf_id.txt
fi
if [ -s matrix_id.txt ]; then
    rm matrix_id.txt
fi

## Indexing of the reference with SAMTools and BWA
## creation of the BED file for the bedcoverage module
## creation of the GATK and picard reference dictionaries

  if [ ! -s "$ref_index" ]; then
	echo -e "Submitting qsub job for BWA reference indexing\n"
    cmd="$BWA index -a is -p ${seq_directory}/${ref} ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N index -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_INDEX_FILE" ]; then
    echo -e "Submitting qsub job for SAMtools reference indexing\n"
    cmd="$SAMTOOLS faidx ${seq_directory}/${ref}.fasta"
    qsub_id=`qsub -N SAM_index -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_DICT" ]; then
    echo -e "Submitting qsub job for ${ref}.dict creation\n"
    cmd="$JAVA $SET_VAR $CREATEDICT R=${seq_directory}/${ref}.fasta O=$REF_DICT"
	qsub_id=`qsub -N PICARD_dict -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
	echo -e "$qsub_id" >> qsub_ids.txt
  fi
  if [ ! -s "$REF_BED" -a qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
	qsub_cat_ids=`cat qsub_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="-hold_jid ${qsub_cat_ids}"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
   echo -e "$qsub_id" >> qsub_ids.txt	
  fi
  if [ ! -s "$REF_BED" -a ! qsub_ids.txt ]; then
    echo -e "Submitting qsub job for BED file construction with BEDTools\n"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    qsub_id=`qsub -N BED_window -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$qsub_id" >> qsub_ids.txt	
  fi

variants_single_sge ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting qsub job for sequence alignment and variant calling for $sequences\n"
        var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N aln_sequences -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
        echo -e "$qsub_array_id" >> qsub_array_ids.txt
	fi
        
fi
if [ ! -s qsub_ids.txt ]; then
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		qsub_array_id=`qsub -N aln_$sequences -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
		echo -e "$qsub_array_id" >> qsub_array_ids.txt
	fi
fi

}

## Variants is the main alignment and variant calling script contained with SPANDx

variants_sge ()
{
if [ -s qsub_ids.txt ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
                var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N aln_${sequences[$i]} -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
                echo -e "$qsub_array_id" >> qsub_array_ids.txt
	        fi
        done
fi
if [ ! -s qsub_ids.txt ]; then
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting qsub job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    	    var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        qsub_array_id=`qsub -N aln_${sequences[$i]} -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
				echo -e "$qsub_array_id" >> qsub_array_ids.txt
	        fi
        done
fi

}


## Matrix is the main comparative genomics section of SPANDx that relies on the outputs from the variant function above
matrix_sge ()			
{	
if [ -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt	
fi
if [ ! -s qsub_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="-hold_jid ${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
		qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
		echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done 
fi
## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes

if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    qsub_cat_ids=`cat clean_vcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
}

## This function will generate a SNP matrix from the Phylo directory assuming all SNP and BAM files have already been linked into these directories
## qsub for indel_merge is currently untested

matrix_final_sge ()			
{
if [ -s qsub_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    qsub_cat_ids=`cat qsub_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
	echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt
fi
if [ ! -s qsub_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting qsub job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Master_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "$qsub_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="-hold_jid ${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes


if [ -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    qsub_cat_ids=`cat mastervcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="-hold_jid ${qsub_cat_ids}"
	echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    echo -e "Submitting qsub job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			qsub_clean_id=`qsub -N clean_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$qsub_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

#####################


if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    qsub_cat_ids=`cat clean_vcf_id.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting qsub job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	qsub_matrix_id=`qsub -N Matrix_vcf -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$qsub_matrix_id" >> matrix_id.txt
fi
}

## This function takes the output of each bedcoverage assessment of the sequence alignment and merges them in a comparative matrix
## This function is run when $strain=all but not when a single strain is analysed i.e. $strain doesn't equal all
merge_BED_sge ()
{
if [ -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    qsub_cat_ids=`cat qsub_array_ids.txt | cut -f3 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="-hold_jid ${qsub_cat_ids}"
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT $depend -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$qsub_BED_id" >> qsub_BED_id.txt
fi
if [ ! -s qsub_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    echo -e "Submitting qsub job for BEDcov merge\n"
    var="ref=$ref,seq_path=$seq_directory"
	qsub_BED_id=`qsub -N BEDcov_merge -j $ERROR_OUT_SGE -m $MAIL -M $ADDRESS -l h_rt=$H_RT -v "$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$qsub_BED_id" >> qsub_BED_id.txt
fi
}

### These variable tests determine which of the above functions need to be run for the SPANDx pipeline

  if [ "$strain" == all ]; then
    variants_sge
	merge_BED_sge
  fi
  if [ "$strain" != all -a "$strain" != none ]; then
    variants_single_sge
  fi
  if [ "$matrix" == yes -a "$strain" != none ]; then
    matrix_sge
  fi
  if [ "$matrix" == yes -a "$strain" == none ]; then
    matrix_final_sge
  fi
fi

################################ SLURM #############################################

if [ "$SCHEDULER" == SLURM ]; then

#clean batch system log files

if [ -s sbatch_ids.txt ]; then
    rm sbatch_ids.txt
fi
if [ -s sbatch_ids2.txt ]; then
    rm sbatch_ids2.txt
fi
if [ -s sbatch_array_ids.txt ]; then
    rm sbatch_array_ids.txt
fi
if [ -s mastervcf_id.txt ]; then
    rm mastervcf_id.txt
fi
if [ -s clean_vcf_id.txt ]; then
	rm clean_vcf_id.txt
fi
if [ -s matrix_id.txt ]; then
    rm matrix_id.txt
fi

## Indexing of the reference with SAMTools and BWA
## creation of the BED file for the bedcoverage module
## creation of the GATK and picard reference dictionaries

  if [ ! -s "$ref_index" ]; then
	echo -e "Submitting sbatch job for BWA reference indexing\n"
    cmd="$BWA index -a is -p ${seq_directory}/${ref} ${seq_directory}/${ref}.fasta"
    sbatch_id=`sbatch --job-name=index --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$sbatch_id" >> sbatch_ids.txt
  fi
  if [ ! -s "$REF_INDEX_FILE" ]; then
    echo -e "Submitting sbatch job for SAMtools reference indexing\n"
    cmd="$SAMTOOLS faidx ${seq_directory}/${ref}.fasta"
    sbatch_id=`sbatch --job-name=SAM_index --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$sbatch_id" >> sbatch_ids.txt
  fi
  if [ ! -s "$REF_DICT" ]; then
    echo -e "Submitting sbatch job for ${ref}.dict creation\n"
    cmd="$JAVA $SET_VAR $CREATEDICT R=${seq_directory}/${ref}.fasta O=$REF_DICT"
	sbatch_id=`sbatch --job-name=PICARD_dict --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
	echo -e "$sbatch_id" >> sbatch_ids.txt
  fi
  if [ ! -s "$REF_BED" -a sbatch_ids.txt ]; then
    echo -e "Submitting sbatch job for BED file construction with BEDTools\n"
	sbatch_cat_ids=`cat sbatch_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="--depend=afterok:${sbatch_cat_ids}"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    sbatch_id=`sbatch --job-name=BED_window --mail-type=$MAIL_SLURM --time=$TIME $depend --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
   echo -e "$sbatch_id" >> sbatch_ids.txt	
  fi
  if [ ! -s "$REF_BED" -a ! sbatch_ids.txt ]; then
    echo -e "Submitting sbatch job for BED file construction with BEDTools\n"
    cmd="$BEDTOOLS makewindows -g $REF_INDEX_FILE -w $window > $REF_BED"
    sbatch_id=`sbatch --job-name=BED_window --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
    echo -e "$sbatch_id" >> sbatch_ids.txt	
  fi

variants_single_slurm ()
{
if [ -s sbatch_ids.txt ]; then
    sbatch_cat_ids=`cat sbatch_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting sbatch job for sequence alignment and variant calling for $sequences\n"
        var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		sbatch_array_id=`sbatch --job-name=aln_sequences --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
        echo -e "$sbatch_array_id" >> sbatch_array_ids.txt
	fi
        
fi
if [ ! -s sbatch_ids.txt ]; then
    if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/$sequences.snps.PASS.vcf ]; then
		echo -e "Submitting sbatch job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    var="seq=$sequences,ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		sbatch_array_id=`sbatch --job-name=aln_$sequences --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
		echo -e "$sbatch_array_id" >> sbatch_array_ids.txt
	fi
fi

}

## Variants is the main alignment and variant calling script contained with SPANDx

variants_slurm ()
{
if [ -s sbatch_ids.txt ]; then
    sbatch_cat_ids=`cat sbatch_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting sbatch job for sequence alignment and variant calling for ${sequences[$i]}\n"
                var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        sbatch_array_id=`sbatch --job-name=aln_${sequences[$i]} --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
                echo -e "$sbatch_array_id" >> sbatch_array_ids.txt
	        fi
        done
fi
if [ ! -s sbatch_ids.txt ]; then
        for (( i=0; i<n; i++ )); do
            if [ ! -s ${PBS_O_WORKDIR}/Outputs/SNPs_indels_PASS/${sequences[$i]}.snps.PASS.vcf ]; then
		        echo -e "Submitting sbatch job for sequence alignment and variant calling for ${sequences[$i]}\n"
	    	    var="seq=${sequences[$i]},ref=$ref,org=$org,strain=$strain,variant_genome=$variant_genome,annotate=$annotate,tech=$tech,pairing=$pairing,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH"
		        sbatch_array_id=`sbatch --job-name=aln_${sequences[$i]} --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/Align_SNP_indel.sh`
				echo -e "$sbatch_array_id" >> sbatch_array_ids.txt
	        fi
        done
fi

}


## Matrix is the main comparative genomics section of SPANDx that relies on the outputs from the variant function above
matrix_slurm ()			
{	
if [ -s sbatch_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    sbatch_cat_ids=`cat sbatch_array_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
    echo -e "Submitting sbatch job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Master_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$sbatch_matrix_id" >> mastervcf_id.txt	
fi
if [ ! -s sbatch_array_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting sbatch job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Master_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/Master_vcf.sh`
	echo -e "$sbatch_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    sbatch_cat_ids=`cat mastervcf_id.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="--depend=afterok:${sbatch_cat_ids}"
	echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
	out=("${sequences_tmp[@]/_1_sequence.fastq.gz/.clean.vcf}")
	bam=("${sequences_tmp[@]/_1_sequence.fastq.gz/.bam}")
	bam_array=("${bam[@]/#/$PBS_O_WORKDIR/Phylo/bams/}")
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${sequences[$i]}.clean.vcf ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
		if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${sequences[$i]}.clean.vcf -a "$indel_merge" == yes ]; then
		cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
		sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
		echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done 
fi
## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes

if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    sbatch_cat_ids=`cat clean_vcf_id.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
    echo -e "Submitting sbatch job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Matrix_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$sbatch_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting sbatch job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Matrix_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$sbatch_matrix_id" >> matrix_id.txt
fi
}

## This function will generate a SNP matrix from the Phylo directory assuming all SNP and BAM files have already been linked into these directories
## sbatch for indel_merge is currently untested

matrix_final_slurm ()			
{
if [ -s sbatch_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    sbatch_cat_ids=`cat sbatch_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
	echo -e "Submitting sbatch job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Master_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "$sbatch_matrix_id" >> mastervcf_id.txt
fi
if [ ! -s sbatch_ids.txt -a ! -s Phylo/out/master.vcf ]; then
    echo -e "Submitting sbatch job for creation of master VCF file\n"
    var="ref=$ref,seq_path=$seq_directory,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Master_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/Master_vcf_final.sh`
	echo -e "$sbatch_matrix_id" >> mastervcf_id.txt
fi

### creates clean.vcf files for SNP calls across all genomes

if [ -s mastervcf_id.txt ]; then
    sbatch_cat_ids=`cat mastervcf_id.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="--depend=afterok:${sbatch_cat_ids}"
	echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt ]; then
    echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/snps/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/out/master.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi

## if indels merge is set to yes creates clean.vcf files for all indels identified across all genomes


if [ -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    sbatch_cat_ids=`cat mastervcf_id.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
	depend="--depend=afterok:${sbatch_cat_ids}"
	echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
	clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi
if [ ! -s mastervcf_id.txt -a "$indel_merge" == yes ]; then
    echo -e "Submitting sbatch job for error checking SNP calls across all genomes\n"
    clean_array=(`find $PBS_O_WORKDIR/Phylo/indels/*.vcf -printf "%f "`)
    out=("${clean_array[@]/.vcf/.clean.vcf}")
    bam_array=(`find $PBS_O_WORKDIR/Phylo/bams/*.bam`)
    n=${#bam_array[@]}
    for (( i=0; i<n; i++ )); do
	    if [ ! -s $PBS_O_WORKDIR/Phylo/out/${out[$i]} ]; then
		    cmd="$JAVA $SET_VAR $GATK -T UnifiedGenotyper -rf BadCigar -R $PBS_O_WORKDIR/${ref}.fasta -I ${bam_array[$i]} -o $PBS_O_WORKDIR/Phylo/indels/out/${out[$i]} -alleles:masterAlleles $PBS_O_WORKDIR/Phylo/indels/out/master_indels.vcf -gt_mode GENOTYPE_GIVEN_ALLELES -out_mode EMIT_ALL_SITES -stand_call_conf 0.0 -glm BOTH -G none"
			sbatch_clean_id=`sbatch --job-name=clean_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export=command="$cmd" "$SCRIPTPATH"/Header.pbs`
			echo -e "$sbatch_clean_id" >> clean_vcf_id.txt
		fi
    done  
fi


if [ -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    sbatch_cat_ids=`cat clean_vcf_id.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
    echo -e "Submitting sbatch job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Matrix_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$sbatch_matrix_id" >> matrix_id.txt
fi
if [ ! -s clean_vcf_id.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    echo -e "Submitting sbatch job for creation of SNP array\n"
    var="ref=$ref,seq_path=$seq_directory,variant_genome=$variant_genome,annotate=$annotate,SCRIPTPATH=$SCRIPTPATH,indel_merge=$indel_merge"
	sbatch_matrix_id=`sbatch --job-name=Matrix_vcf --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/SNP_matrix.sh`
	echo -e "$sbatch_matrix_id" >> matrix_id.txt
fi
}

## This function takes the output of each bedcoverage assessment of the sequence alignment and merges them in a comparative matrix
## This function is run when $strain=all but not when a single strain is analysed i.e. $strain doesn't equal all
merge_BED_slurm ()
{
if [ -s sbatch_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    sbatch_cat_ids=`cat sbatch_array_ids.txt | cut -f4 -d ' ' | sed -e 's/$/,/' | tr -d '\n' | sed -e 's/,$//'`
    depend="--depend=afterok:${sbatch_cat_ids}"
    echo -e "Submitting sbatch job for BEDcov merge\n"
    var="seq_path=$seq_directory"
	sbatch_BED_id=`sbatch --job-name=BEDcov_merge --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME $depend --export="$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$sbatch_BED_id" >> sbatch_BED_id.txt
fi
if [ ! -s sbatch_array_ids.txt -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Bedcov_merge.txt ]; then
    echo -e "Submitting sbatch job for BEDcov merge\n"
    var="ref=$ref,seq_path=$seq_directory"
	sbatch_BED_id=`sbatch --job-name=BEDcov_merge --mem=$SLURM_MEM --mail-type=$MAIL_SLURM --time=$TIME --export="$var" "$SCRIPTPATH"/BedCov_merge.sh`
	#echo -e "$sbatch_BED_id" >> sbatch_BED_id.txt
fi
}

### These variable tests determine which of the above functions need to be run for the SPANDx pipeline

  if [ "$strain" == all ]; then
    variants_slurm
	merge_BED_slurm
  fi
  if [ "$strain" != all -a "$strain" != none ]; then
    variants_single_slurm
  fi
  if [ "$matrix" == yes -a "$strain" != none ]; then
    matrix_slurm
  fi
  if [ "$matrix" == yes -a "$strain" == none ]; then
    matrix_final_slurm
  fi
fi

exit 0
Status API Training Shop Blog About
© 2015 GitHub, Inc. Terms Privacy Security Contact
