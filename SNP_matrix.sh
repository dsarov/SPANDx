#!/bin/bash
#$ -S /bin/bash
#PBS -S /bin/sh
#PBS -o logs/SNP_matrix.txt
#$ -cwd

#########################################################################
# The following script will combine all vcf files in your GATK directory into a master VCF file
# Combined snp calls from the master.vcf are interrogated across all bam files at these locations to verify snp identity
# Finally the clean vcf files are concatenated and converted into a matrix for phylogeny programs
# Version 1.5
# Written by D. Sarovich
# derek.sarovich@menzies.edu.au
#
# Version history
# Version 1.1-1.2 
# Added post processing of SNP matrix to remove lines containing 1 and -1 and combined individual samples and IDs into a single matrix
# Some additional file cleanup 
# V1.2-1.3
# annotated indels and snps are now tab delimited. Tri-allelic calls should be handled better in the indels. A previous bug was causing the nucleotide output to include 2/2
# V1.4
# added SGE header and -cwd for correct processing of log files
# V1.5
# added version control to snpEff and fixed processing of tri-allelic snps and indels in the annotated snps/indels file
#
#########################################################################

#source variables
source "$SCRIPTPATH"/SPANDx.config

if [ ! $PBS_O_WORKDIR ]
    then
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
#bgzip and tabix clean vcf files. Combines clean vcf files into an array
#########################################################################

if [ ! -s $PBS_O_WORKDIR/Phylo/out/out.vcf ]; then
    log_eval $PBS_O_WORKDIR/Phylo/out "for f in *.clean.vcf; do $BGZIP -c \$f > \${f}.gz; done;"
    log_eval $PBS_O_WORKDIR/Phylo/out "for f in *.clean.vcf.gz; do $TABIX -p vcf \$f; done;"
    vcf_array=`(find $PBS_O_WORKDIR/Phylo/out/*.vcf.gz -printf "%f ")`
    log_eval $PBS_O_WORKDIR/Phylo/out "$VCFMERGE ${vcf_array[*]} >out.vcf"   # This command results in a broken pipe if the number of genomes doesn't match the number of SNP/VCF files
  
  ## to do
  ## include some verification process to make sure the arrays match for matrix creation and error output if this isn't the case
  
    log_eval $PBS_O_WORKDIR/Phylo/out "$VCFTOOLS --vcf out.vcf --out merge --012"
  else 
    echo -e "out.vcf has already been created\n\n"
fi
if [ -s $PBS_O_WORKDIR/Phylo/out/merge.012.indv -a ! -s $PBS_O_WORKDIR/Phylo/out/merge.012.indv.trans ] 
  then
  
    cd $PBS_O_WORKDIR/Phylo/out
    awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' merge.012 >merge.012.trans
    echo -e ${ref} > ref.snp
	cat ref.snp merge.012.indv > merge.012.indv.ref
	awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' merge.012.indv.ref >merge.012.indv.trans
	cd $PBS_O_WORKDIR
fi

#########################################################################
#bgzip and tabix clean vcf files. Combines clean vcf files into an array if indels merge is set to yes
#########################################################################  
if [ ! -s $PBS_O_WORKDIR/Phylo/indels/out/out.vcf -a "$indel_merge" == yes ]; then
    log_eval $PBS_O_WORKDIR/Phylo/indels/out "for f in *.clean.vcf; do $BGZIP -c \$f > \${f}.gz; done;" 
	log_eval $PBS_O_WORKDIR/Phylo/indels/out "for f in *.clean.vcf.gz; do $TABIX -p vcf \$f; done;"
	vcf_array=`(find $PBS_O_WORKDIR/Phylo/indels/out/*.vcf.gz -printf "%f ")`
    log_eval $PBS_O_WORKDIR/Phylo/indels/out "$VCFMERGE ${vcf_array[*]} >out.vcf"
	log_eval $PBS_O_WORKDIR/Phylo/indels/out "$VCFTOOLS --vcf out.vcf --out merge --012"
	else 
    echo -e "indels out.vcf has already been created\n\n"
fi
if [ -s $PBS_O_WORKDIR/Phylo/indels/out/merge.012.indv -a ! -s $PBS_O_WORKDIR/Phylo/indels/out/merge.012.indv.trans ] 
  then
    cd $PBS_O_WORKDIR/Phylo/indels/out
    awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' merge.012 >merge.012.trans
	echo -e ${ref} > ref.indel
	cat ref.indel merge.012.indv > merge.012.indv.ref
    awk 'BEGIN {FS=OFS=" "}{for (i=1;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=1;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' merge.012.indv.ref >merge.012.indv.trans
	cd $PBS_O_WORKDIR
fi
###########################################################################
# Merge the individual files created with VCF tools, remove the lines containing 1 and -1 and clean-up files
############################################################################

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    cd $PBS_O_WORKDIR/Phylo/out/
# merge the chromosomes and position so they are in the chromosomes_position format
    sed 's/\t/_/' merge.012.pos > merge.012.pos_merged
# remove the header from the out.vcf file and the header from the 012 matrix
    grep -v '#' out.vcf > out.vcf.headerless
    tail -n +2 merge.012.trans > merge.012.trans.top_rem
# print the vcf.out positions and merge to the same format as merge.012.pos_merged
    awk '{print $1, $2}' out.vcf.headerless > out.vcf.pos
    sed 's/ /_/' out.vcf.pos > out.vcf.pos_merged
# print the reference and alternate allele columns from the merged VCF file	
    awk '{print $4, $5}' out.vcf.headerless > out.vcf.ref_alt_AGCT
#paste in positions
    paste out.vcf.pos_merged out.vcf.ref_alt_AGCT > vcf.pos.alleles.AGCT
## filter vcf positions for those present in the final 012 matrix and output only those lines
    awk 'FNR==NR{ a[$1]=$0;next } ($1 in a)' merge.012.pos_merged vcf.pos.alleles.AGCT > filtered.pos.AGCT
## merge the 012 matrix
    paste filtered.pos.AGCT merge.012.trans.top_rem > pos.alleles.AGCT.012 
#remove 1 and -1
    grep -v ' 1' <pos.alleles.AGCT.012 | grep -v ' -1' | sed '/\t-1/d' | sed '/\t1/d' >t1
    awk '  { for (i=4; i<=NF; i++) {if ($i == 0) $i=$2; if ($i == 2) $i=$3 }}; {print $0} ' t1 > t2 
    cut -d " " -f-2,4- t2>t3
	x=`cat t3 | wc -l`
    y=`cat merge.012.indv.ref | wc -l`
    taxa=`cat merge.012.indv.trans`
    grid=`cat t3`
   #nexus creation for PAUP
    echo -e "\n#nexus\nbegin data;\ndimensions ntax=$y nchar=$x;\nformat symbols=\"AGCT\" gap=. transpose;\ntaxlabels $taxa;\nmatrix\n$grid\n;\nend;" > Ortho_SNP_matrix.nex
	mv $PBS_O_WORKDIR/Phylo/out/Ortho_SNP_matrix.nex $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex
fi
if [ ! -s $PBS_O_WORKDIR/Phylo/out/Phylo_RAxML.nex -a ! -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex ]; then 
	awk 'BEGIN {FS=" "; OFS=""}{for (i=2;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=2;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' t3 > t4
	paste merge.012.indv.ref t4 > t5
	x=`cat t3 | wc -l` #nchar
	y=`cat merge.012.indv.ref | wc -l` #ntaxa
	taxa_and_grid=`cat t5`
	echo -e "$y $x\n$taxa_and_grid" > Ortho_SNP_matrix_RAxML.nex
	mv $PBS_O_WORKDIR/Phylo/out/Ortho_SNP_matrix_RAxML.nex $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_RAxML.nex
fi	

############################
## construct a SNP matrix including tri and tetra allelic SNPs

if [ $tri_tetra_allelic == yes ]; then
    cd $PBS_O_WORKDIR/Phylo/out/
	grep -v '#' out.vcf | grep -v '\./\.' | grep -v '0/1' | grep -v '0/2' | grep -v '0/3'| grep -v '1/2' | grep -v '1/3'| grep -v '2/3' >headerless_vcf

	
	#replace genotype calls with nucleotide information
	
   awk '  { for (i=10; i<=NF; i++) 
           { if ($i ~ /0\/0/) $i=$4; 
            if ($i ~ /1\/1/) $i=substr($5, 1, 1);
            if ($i ~ /2\/2/) $i=substr($5, 3, 1);
            if ($i ~ /3\/3/) $i=substr($5, 5, 1);
             }}
            {print $0} ' headerless_vcf > clean.out.tri_allelic.vcf


    grep '#CHROM' out.vcf > simple.header
	
	cat clean.out.tri_allelic.vcf | awk '{print $4}' > ref_column
	
	#PAUP nex file
	
	awk '{print $1,$2}' clean.out.tri_allelic.vcf | sed 's/ /_/' > locations
	cat clean.out.tri_allelic.vcf | cut -d' ' -f10- | sed 's/ /\t/g' >calls
	paste locations ref_column calls > matrix
	grid=`cat matrix`
    x=`cat clean.out.tri_allelic.vcf | wc -l` #nchars
    y=`cat merge.012.indv.ref | wc -l` #ntaxa
    taxa=`cat merge.012.indv.trans`
    
    echo -e "\n#nexus\nbegin data;\ndimensions ntax=$y nchar=$x;\nformat symbols=\"AGCT\" gap=. transpose;\ntaxlabels $taxa;\nmatrix\n$grid\n;\nend;" > Ortho_SNP_matrix_poly_allelic.nex
    mv $PBS_O_WORKDIR/Phylo/out/Ortho_SNP_matrix_poly_allelic.nex $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_poly_allelic.nex
	
	#RAXML
	
	awk 'BEGIN {FS=" "; OFS=""}{for (i=2;i<=NF;i++){arr[NR,i]=$i; if(big <= NF) big=NF;}}END {for(i=2;i<=big;i++){for(j=1;j<=NR;j++){printf("%s%s",arr[j,i],(j==NR?"":OFS));} print "";}}' matrix > matrix_RAXML
	paste merge.012.indv.ref matrix_RAXML > t5
	taxa_and_grid=`cat t5`
	echo -e "$y $x\n$taxa_and_grid" > Ortho_SNP_matrix_poly_allelic_RAxML.nex
	mv $PBS_O_WORKDIR/Phylo/out/Ortho_SNP_matrix_poly_allelic_RAxML.nex $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix_poly_allelic_RAxML.nex
	
	rm matrix matrix_RAXML
	
	
fi


##################indels ##############
if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex -a "$indel_merge" == yes ]; then
    cd $PBS_O_WORKDIR/Phylo/indels/out
    sed 's/\t/_/' merge.012.pos > merge.012.pos_merged
    grep -v '#' out.vcf > out.vcf.headerless
    tail -n +2 merge.012.trans > merge.012.trans.top_rem
    awk '{print $1, $2}' out.vcf.headerless > out.vcf.pos
    sed 's/ /_/' out.vcf.pos > out.vcf.pos_merged
    awk '{print $4, $5}' out.vcf.headerless > out.vcf.ref_alt_AGCT
    paste out.vcf.pos_merged out.vcf.ref_alt_AGCT > vcf.pos.alleles.AGCT
    awk 'FNR==NR{ a[$1]=$0;next } ($1 in a)' merge.012.pos_merged vcf.pos.alleles.AGCT > filtered.pos.AGCT
    paste filtered.pos.AGCT merge.012.trans.top_rem > pos.alleles.AGCT.012 
    grep -v ' 1' <pos.alleles.AGCT.012 | grep -v ' -1' | sed '/\t-1/d' | sed '/\t1/d' >t1
    awk '  { for (i=4; i<=NF; i++) {if ($i == 0) $i=$2; if ($i == 2) $i=$3 }}; {print $0} ' t1 > t2 
    cut -d " " -f-2,4- t2>t3
    x=`cat t3 | wc -l`
    y=`cat merge.012.indv | wc -l`
    taxa=`cat merge.012.indv.trans`
    grid=`cat t3`
   #nexus creation for PAUP
    echo -e "\n#nexus\nbegin data;\ndimensions ntax=$y nchar=$x;\nformat symbols=\"AGCT\" gap=. transpose;\ntaxlabels $taxa;\nmatrix\n$grid\n;\nend;" > indel_matrix.nex
    mv $PBS_O_WORKDIR/Phylo/indels/out/indel_matrix.nex $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex
fi


###############################################
##
## These steps will take the merged SNP outputs, annotate them and reformat the data into a tab delimited txt file for importation into excel
## output from these steps should be a tab delimited file
##
## Annotation of SNPs using snpEff
##
##################################################

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs_annotated.txt -a "$annotate" == yes ]; then
    if [ ! -d $PBS_O_WORKDIR/Phylo/annotated ]; then
	    mkdir $PBS_O_WORKDIR/Phylo/annotated
	fi
	cp $PBS_O_WORKDIR/Phylo/out/out.vcf $PBS_O_WORKDIR/Phylo/annotated/out.vcf
	cp $PBS_O_WORKDIR/Phylo/out/merge.012.indv.trans $PBS_O_WORKDIR/Phylo/annotated/merge.012.indv.trans
	cd $PBS_O_WORKDIR/Phylo/annotated
	
	#SnpEff version control. Versions post 4.1 use a slightly different format for variants
	
	SnpEff_version=`($JAVA $SET_VAR $SNPEFF -h 2>&1 | head -n1 |awk '{ print $4 }' | awk '{print substr($1,0,3)}' | bc)`
	if [ "$(echo $SnpEff_version '>=' 4.1 | bc -l)" -eq 1 ]; then
		
		echo -e "Version of snpEff is 4.1 or greater\n"
		echo -e "Version of SnpEff is $SnpEff_version"
		log_eval $PBS_O_WORKDIR/Phylo/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -formatEff -v $variant_genome $PBS_O_WORKDIR/Phylo/annotated/out.vcf > $PBS_O_WORKDIR/Phylo/annotated/out.annotated.vcf"
	else
		echo -e "Version of snpEff is less than 4.1\n"
		log_eval $PBS_O_WORKDIR/Phylo/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -v $variant_genome $PBS_O_WORKDIR/Phylo/annotated/out.vcf > $PBS_O_WORKDIR/Phylo/annotated/out.annotated.vcf"
	fi	
		
	rm snpEff_genes.txt
	rm snpEff_summary.html
		
	#Parse annotated vcf and merged vcf into multiple files for processing
	#effects.mrg includes all annotations
	#out.pos.mrg includes every SNP position
	#bases.raw contains the nucleotide information for the mutations
	#--------------

	#remove headers from annotated vcf and out.vcf
	grep -v '#' out.annotated.vcf > out.annotated.vcf.headerless
	grep -v '#' out.vcf > out.vcf.headerless

	### effects.mrg ####
	#print effects into separate file and process for merge 
	awk '{print $8 }' out.annotated.vcf.headerless > effects.temp
	awk ' BEGIN {FS=";" } ; {print $6}' effects.temp > effects
	sed -i 's/EFF=//' effects
	sed -i 's/(/ /g' effects
	sed -i 's/|/ /g' effects
	sed -i 's/UPSTREAM MODIFIER /UPSTREAM MODIFIER - /g' effects
	cut -d " " -f -8 effects > effects.mrg
	rm effects effects.temp

	### out.pos.mrg ###	
	awk '{print $1, $2}' out.annotated.vcf.headerless > out.pos.mrg
	#merge the chromosome name and chromosome location into a single field
	sed -i 's/ /_/' out.pos.mrg

	
	### bases.raw ####
	##parse out base information
	awk '{print $4, $5 }' out.vcf.headerless > bases.raw
	sed -i 's/,/ /g' bases.raw
	sed -i 's/$/ N N/' bases.raw
	cut -d " " -f -4 bases.raw > bases.raw.parsed
	sed -i 's/ /\t/g' bases.raw.parsed
		
		
	#-------------------------------------
	# print just the SNP/genomic information column 10 until the end
	cut -f 10- out.vcf.headerless > out.snps.raw

	awk '{for (i=1; i<=NF; i++) { $i=(substr($i, 1, 3))}; print $0 }' out.snps.raw > out.snps.clean
	##remove blank lines at end of out.snps.clean

	#replace SNP information i.e. 1/1, 1/2 etc with 0123 matrix
	
	awk ' { for (i=1; i<=NF; i++) 
		{if ($i == "1/1") $i=2; 
		 if ($i == "./.") $i="."; 
		 if ($i == "0/0") $i=0; 
		 if ($i == "0/1") $i="?"; 
		 if ($i == "2/2") $i="3"; 
		 if ($i == "1/2") $i="?"; 
		 if ($i == "0/2") $i="?";
		 if ($i == "2/1") $i="?";
		 if ($i == "3/3") $i="4";
		 if ($i == "0/3") $i="?";
		 if ($i == "1/3") $i="?";
		 if ($i == "2/3") $i="?";
		 if ($i == "3/2") $i="?";	
		 if ($i == "3/1") $i="?"}}; 
		{print $0} ' out.snps.clean > out.matrix

	#add the reference SNP call to each line and remove spaces to create the 012.mrg for the final matrix
	sed -i -e 's/^/0 /' out.matrix
	sed 's/ //g' out.matrix > 012.mrg

	##convert to tab separated
	sed -i 's/ /\t/g' out.matrix
    # merge base calls and 012 matrix 
	
    paste -d "\t" bases.raw.parsed out.matrix > bases.tab
	
	#replace 012 matrix with base calls
	awk '  { for (i=4; i<=NF; i++) 
	       {if ($i == 0) $i=$1; 
		    if ($i == 2) $i=$2;
            if ($i == 3) $i=$3;			
			if ($i == 4) $i=$4 }}; {print $0} ' bases.tab > bases.tab.tmp
	
	## removes temp columns for final base call file
	cut -d " " -f 5- bases.tab.tmp > bases.mrg

	## paste all temp .mrg files into a headerless effects file
	paste -d " " out.pos.mrg bases.mrg 012.mrg effects.mrg > headerless.effects.mrg

	#creation of the header and merging of files

	echo -e "Location" >> location.mrg
	echo -e "Binary_code Effect Impact Functional_Class Codon_change Amino_Acid_change Amino_Acid_Length Gene_name Biotype" >> header.mrg
	paste -d " " location.mrg merge.012.indv.trans header.mrg > final.header.mrg
	cat final.header.mrg headerless.effects.mrg > All_SNPs_annotated.txt
	sed -i 's/ /\t/g' All_SNPs_annotated.txt
	
	cleanup () 
	{
	rm out.pos.mrg bases.mrg 012.mrg effects.mrg
	rm out.snps.raw out.snps.clean out.matrix bases.raw
	rm out.vcf.headerless out.annotated.vcf.headerless
	rm bases.tab.tmp #check bases.mrg.tmp here
	rm final.header.mrg header.mrg location.mrg headerless.effects.mrg
	rm out.vcf merge.012.indv.trans
	rm out.annotated.vcf
	}

	##run cleanup
	cleanup
	mv $PBS_O_WORKDIR/Phylo/annotated/All_SNPs_annotated.txt $PBS_O_WORKDIR/Outputs/Comparative/All_SNPs_annotated.txt		
fi

########################################
##                                    ##
## Annotation of indels using snpEff  ##
##                                    ##
########################################

if [ ! -s $PBS_O_WORKDIR/Outputs/Comparative/All_indels_annotated.txt -a "$annotate" == yes -a "$indel_merge" == yes ]; then
    if [ ! -d $PBS_O_WORKDIR/Phylo/annotated ]; then
	    mkdir $PBS_O_WORKDIR/Phylo/annotated
	fi
	cp $PBS_O_WORKDIR/Phylo/indels/out/out.vcf $PBS_O_WORKDIR/Phylo/annotated/out_indels.vcf
	cd $PBS_O_WORKDIR/Phylo/annotated 
	
	
	#SnpEff version control. Versions post 4.1 use a slightly different format for variants
	
	
	

	
	SnpEff_version=`($JAVA $SET_VAR $SNPEFF -h 2>&1 | head -n1 |awk '{ print $4 }' | awk '{print substr($1,0,3)}' | bc)`
	if [ "$(echo $SnpEff_version '>=' 4.1 | bc -l)" -eq 1 ]; then
		
		echo -e "Version of snpEff is 4.1 or greater\n"
		echo -e "Version of SnpEff is $SnpEff_version"
		log_eval $PBS_O_WORKDIR/Phylo/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -formatEff -v $variant_genome $PBS_O_WORKDIR/Phylo/annotated/out_indels.vcf > $PBS_O_WORKDIR/Phylo/annotated/out_indels.annotated.vcf"
	else
		echo -e "Version of snpEff is less than 4.1\n"
		log_eval $PBS_O_WORKDIR/Phylo/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -v $variant_genome $PBS_O_WORKDIR/Phylo/annotated/out_indels.vcf > $PBS_O_WORKDIR/Phylo/annotated/out_indels.annotated.vcf"
	fi	
    
	rm snpEff_genes.txt
	rm snpEff_summary.html
	cp $PBS_O_WORKDIR/Phylo/indels/out/merge.012.indv.trans $PBS_O_WORKDIR/Phylo/annotated/merge.012.indv.trans ## This file should be removed in the cleanup section
  
    #Parse annotated vcf and merged vcf into multiple files for processing
    #effects.mrg includes all annotations
    #out.pos.mrg includes every indel position
    #bases.raw contains the nucleotide information for the mutations
    #--------------

    #remove headers from annotated vcf and out.vcf
    grep -v '#' out_indels.annotated.vcf > out.annotated.vcf.headerless
    grep -v '#' out_indels.vcf > out.vcf.headerless

    ### effects.mrg ####
    #print effects into separate file and process for merge 
    awk '{print $8 }' out.annotated.vcf.headerless > effects.temp
    awk ' BEGIN {FS=";" } ; {print $6}' effects.temp > effects
    sed -i 's/EFF=//' effects
    sed -i 's/(/|/g' effects
    sed -i 's/|/\t/g' effects
	sed -i 's/)//g' effects
	sed 's/,/\t/g' effects > effects.mrg
   # sed -i 's/UPSTREAM MODIFIER /UPSTREAM MODIFIER - /g' effects
  #  cut -d " " -f -8 effects > effects.mrg
    rm effects effects.temp

    ### out.pos.mrg ###	
    awk '{print $1, $2}' out.annotated.vcf.headerless > out.pos.mrg
    #merge the chromosome name and chromosome location into a single field
    sed -i 's/ /_/' out.pos.mrg

	
    ### bases.raw ####
    ##parse out base information
    awk ' {print $4, $5 }' out.vcf.headerless > bases.raw
	sed -i 's/,/ /g' bases.raw
	sed -i 's/$/ N N N N N N N/' bases.raw
	cut -d " " -f -9 bases.raw > bases.raw.parsed
	sed -i 's/ /\t/g' bases.raw.parsed

    #-------------------------------------
    # print just the indel information column 10 until the end
    cut -f 10- out.vcf.headerless > out.indels.raw

	#cleanup the raw indel output to just contain the 1/1, 1/2 etc vcf file calls
    awk ' {for (i=1; i<=NF; i++) { $i=(substr($i, 1, 3))}; print $0 }' out.indels.raw > out.indels.clean
    

    #replace indel information i.e. 1/1, 1/2 etc with 0123 matrix and convert to tab separated
	#This section will still have issues with triallelic calls
	
    awk ' { for (i=1; i<=NF; i++) 
	      {if ($i == "1/1") $i="2"; 
		  else if ($i == "./.") $i="."; 
		  else if ($i == "0/0") $i="0"; 
		  else if ($i == "0/1") $i="?"; 
		  else if ($i == "2/2") $i="3"; 
		  else if ($i == "1/2") $i="?"; 
		  else if ($i == "2/1") $i="?"; 
		  else if ($i == "0/2") $i="?";
		  else if ($i == "3/3") $i="4";
		  else if ($i == "0/3") $i="?";
		  else if ($i == "1/3") $i="?";
		  else if ($i == "2/3") $i="?";
		  else if ($i == "3/2") $i="?";	
		  else if ($i == "3/1") $i="?";
		  else if ($i == "4/4") $i="5";
		  else if ($i == "5/5") $i="6";
		  else if ($i == "6/6") $i="7";
		  else if ($i == "7/7") $i="8";
		  else if ($i == "8/8") $i="9";
		  else $i="?"}}; {print $0} ' out.indels.clean > out.matrix
	
	
	
    #add the reference SNP call to each line and remove spaces to create the 012.mrg for the final matrix
    sed -i -e 's/^/0 /' out.matrix
    sed 's/ //g' out.matrix > 012.mrg
    
	
	##convert to tab separated
	sed -i 's/ /\t/g' out.matrix
    # merge base calls and 012 matrix 
    paste -d "\t" bases.raw.parsed out.matrix > bases.tab
    #replace 012 matrix with base calls
    awk ' { for (i=9; i<=NF; i++) 
	{if ($i == 0) $i=$1; 
	 if ($i == 2) $i=$2; 
	 if ($i == 3) $i=$3;
	 if ($i == 4) $i=$4;
	 if ($i == 5) $i=$5;
	 if ($i == 6) $i=$6;
	 if ($i == 7) $i=$7;
	 if ($i == 8) $i=$8;
	 if ($i == 9) $i=$9}}; {print $0} ' bases.tab > bases.mrg.tmp
	 
	 
    ## removes temp columns for final base call file
    cut -d " " -f 10- bases.mrg.tmp > bases.mrg

    ## paste all temp .mrg files into a headerless effects file
    paste -d " " out.pos.mrg bases.mrg 012.mrg effects.mrg > headerless.effects.mrg

    #creation of the header and merging of files

    echo -e "Location" >> location.mrg
    echo -e "Binary_code Effect Impact Functional_Class Codon_change Amino_Acid_change Amino_Acid_Length Gene_Name Transcript_Biotype Gene_Coding Transcript_ID Exon_Rank Genotype_nember Extra_effects_Warnings" >> header.mrg
    paste -d " " location.mrg merge.012.indv.trans header.mrg > final.header.mrg
    cat final.header.mrg headerless.effects.mrg > All_indels_annotated.txt
	sed -i 's/ /\t/g' All_indels_annotated.txt
	
    cleanup () 
    {
    rm out.pos.mrg bases.mrg 012.mrg effects.mrg
    rm out.indels.raw out.indels.clean out.matrix bases.raw
    rm out.vcf.headerless out.annotated.vcf.headerless
    rm bases.tab bases.raw.parsed bases.mrg.tmp merge.012.indv.trans
    rm final.header.mrg header.mrg location.mrg headerless.effects.mrg
	rm out_indels.annotated.vcf out_indels.vcf
    }

    ##run cleanup
    cleanup

	mv $PBS_O_WORKDIR/Phylo/annotated/All_indels_annotated.txt $PBS_O_WORKDIR/Outputs/Comparative/All_indels_annotated.txt
fi

##file clean-up from above matrix steps
if [ -s $PBS_O_WORKDIR/Outputs/Comparative/Ortho_SNP_matrix.nex ]; then
    [ -f $PBS_O_WORKDIR/Phylo/out/t1 ] && rm $PBS_O_WORKDIR/Phylo/out/t1
    [ -f $PBS_O_WORKDIR/Phylo/out/t2 ] && rm $PBS_O_WORKDIR/Phylo/out/t2
    [ -f $PBS_O_WORKDIR/Phylo/out/t3 ] && rm $PBS_O_WORKDIR/Phylo/out/t3
    [ -f $PBS_O_WORKDIR/Phylo/out/t4 ] && rm $PBS_O_WORKDIR/Phylo/out/t4
    [ -f $PBS_O_WORKDIR/Phylo/out/t5 ] && rm $PBS_O_WORKDIR/Phylo/out/t5
    [ -f $PBS_O_WORKDIR/Phylo/out/pos.alleles.AGCT.012 ] && rm $PBS_O_WORKDIR/Phylo/out/pos.alleles.AGCT.012
	[ -f $PBS_O_WORKDIR/Phylo/out/filtered.pos.AGCT ] && rm $PBS_O_WORKDIR/Phylo/out/filtered.pos.AGCT
	[ -f $PBS_O_WORKDIR/Phylo/out/vcf.pos.alleles.AGCT ] && rm $PBS_O_WORKDIR/Phylo/out/vcf.pos.alleles.AGCT
    [ -f $PBS_O_WORKDIR/Phylo/out/out.vcf.ref_alt_AGCT ] && rm $PBS_O_WORKDIR/Phylo/out/out.vcf.ref_alt_AGCT
    [ -f $PBS_O_WORKDIR/Phylo/out/out.vcf.pos_merged ] && rm $PBS_O_WORKDIR/Phylo/out/out.vcf.pos_merged
	[ -f $PBS_O_WORKDIR/Phylo/out/out.vcf.pos ] && rm $PBS_O_WORKDIR/Phylo/out/out.vcf.pos
	[ -f $PBS_O_WORKDIR/Phylo/out/out.vcf.headerless ] && rm $PBS_O_WORKDIR/Phylo/out/out.vcf.headerless
	[ -f $PBS_O_WORKDIR/Phylo/out/merge.012.trans.top_rem ] && rm $PBS_O_WORKDIR/Phylo/out/merge.012.trans.top_rem
	[ -f $PBS_O_WORKDIR/Phylo/out/merge.012.pos_merged ] && rm $PBS_O_WORKDIR/Phylo/out/merge.012.pos_merged
fi


##file clean-up from above matrix steps for indels
if [ -s $PBS_O_WORKDIR/Outputs/Comparative/indel_matrix.nex ]; then
    [ -f $PBS_O_WORKDIR/Phylo/indels/out/t1 ] && rm $PBS_O_WORKDIR/Phylo/indels/out/t1
    [ -f $PBS_O_WORKDIR/Phylo/indels/out/t2 ] && rm $PBS_O_WORKDIR/Phylo/indels/out/t2
    [ -f $PBS_O_WORKDIR/Phylo/indels/out/t3 ] && rm $PBS_O_WORKDIR/Phylo/indels/out/t3
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/pos.alleles.AGCT.012 ] && rm $PBS_O_WORKDIR/Phylo/indels/out/pos.alleles.AGCT.012
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/filtered.pos.AGCT ] && rm $PBS_O_WORKDIR/Phylo/indels/out/filtered.pos.AGCT
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/vcf.pos.alleles.AGCT ] && rm $PBS_O_WORKDIR/Phylo/indels/out/vcf.pos.alleles.AGCT
    [ -f $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.ref_alt_AGCT ] && rm $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.ref_alt_AGCT
    [ -f $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.pos_merged ] && rm $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.pos_merged
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.pos ] && rm $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.pos
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.headerless ] && rm $PBS_O_WORKDIR/Phylo/indels/out/out.vcf.headerless
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/merge.012.trans.top_rem ] && rm $PBS_O_WORKDIR/Phylo/indels/out/merge.012.trans.top_rem
	[ -f $PBS_O_WORKDIR/Phylo/indels/out/merge.012.pos_merged ] && rm $PBS_O_WORKDIR/Phylo/indels/out/merge.012.pos_merged
fi


echo "SPANDx has finished"

sleep 20
exit 0
