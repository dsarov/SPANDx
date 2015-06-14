#!/bin/bash
#PBS -S /bin/sh
#$ -S /bin/bash
#$ -cwd



############################################################################
# 
# Thank you for using SPANDx. SPANDx is meant to be run using the main script SPANDx.sh
# For help please refer to the SPANDx manual or the main script
# If you want to customise the programs run in the pipeline please refer to the SPANDx 
# manual and the comments included in each stage of the program in this script
# Almost all of the commands below test for completion of the next stage in the pipeline before the commands are executed
# Make sure the file checks make sense if customising different sections of the pipeline
#
# Written by Derek Sarovich - Menzies School of Health Research
# Version 2.1
# v2.0-2.1 added SGE header
#
# 
#
############################################################################
#### Variables required seq ref org strain seq_path annotate pairing tech window ####

#source variables
source "$SCRIPTPATH"/SPANDx.config
source "$SCRIPTPATH"/GATK.config

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

###file variable list

REF_FILE=$seq_path/${ref}
REF_FASTA=$seq_path/${ref}.fasta
READ1_FILE=$seq_path/${seq}_1_sequence.fastq.gz
READ2_FILE=$seq_path/${seq}_2_sequence.fastq.gz
SAI1=$seq_path/${seq}/${seq}_1.sai
SAI2=$seq_path/${seq}/${seq}_2.sai
SAM=$seq_path/${seq}/unique/${seq}.sam
BAM=$seq_path/${seq}/unique/${seq}.tmp
TMP_BAM1=$seq_path/${seq}/unique/${seq}_tmp1.bam
TMP_BAM2=$seq_path/${seq}/unique/${seq}_tmp2.bam
TMP_BAM3=$seq_path/${seq}/unique/${seq}_tmp3.bam
UNMAPPED=$seq_path/${seq}/unique/${seq}.unmapped
BAM_INDEX_FILE=$seq_path/${seq}/unique/${seq}.bam.bai
BAM_UNIQUE_FILE=$seq_path/${seq}/unique/${seq}
BAM_UNIQUE_INDEX_FILE=$seq_path/${seq}/unique/${seq}.bam.bai
BAM_DUPLICATESREM_FILE=$seq_path/${seq}/unique/${seq}.dups.bam
BAM_DUPLICATESREM_FILE_INDEX=$seq_path/${seq}/unique/${seq}.dups.bam.bai
DUP_METRICS=$seq_path/${seq}/unique/dup_metrics.txt
GATK_DEPTH_FILE=$seq_path/${seq}/unique/${seq}.sample_summary
GATK_TARGET_FILE=$seq_path/${seq}/unique/${seq}.bam.list
GATK_REALIGNED_BAM=$seq_path/${seq}/unique/${seq}.realigned.bam
GATK_RAW_SNPS_FILE=$seq_path/${seq}/unique/${seq}.snps.raw.vcf
GATK_RAW_INDELS_FILE=$seq_path/${seq}/unique/${seq}.indels.raw.vcf
GATK_FILTER_SNPS_FILE=$seq_path/${seq}/unique/${seq}.snps.filter.vcf
GATK_FILTER_INDELS_FILE=$seq_path/${seq}/unique/${seq}.indels.filter.vcf
GATK_PASS_SNPS_FILE=$seq_path/${seq}/unique/${seq}.snps.PASS.vcf
GATK_PASS_INDELS_FILE=$seq_path/${seq}/unique/${seq}.indels.PASS.vcf
ANNOTATED_SNPS=$seq_path/${seq}/unique/annotated/${seq}.snps.PASS.vcf.annotated
ANNOTATED_INDELS=$seq_path/${seq}/unique/annotated/${seq}.indels.PASS.vcf.annotated
GATK_SNP_HEAD=$seq_path/${seq}/unique/${seq}.snp.head.vcf
GATK_INDEL_HEAD=$seq_path/${seq}/unique/${seq}.indel.head.vcf
BED_WINDOW_FILE=$REF_FILE.bed
BED_COV_FILE=$seq_path/${seq}/unique/${seq}.bedcov
GATK_SNP_AMB=$seq_path/${seq}/unique/${seq}.snps.AMB.vcf
GATK_INDELS_AMB=$seq_path/${seq}/unique/${seq}.indels.AMB.vcf
GATK_SNP_HEAD2=$seq_path/${seq}/unique/${seq}.snp.head2.vcf
GATK_SNP_FAIL=$seq_path/${seq}/unique/${seq}.snps.FAIL.vcf
GATK_INDELS_HEAD2=$seq_path/${seq}/unique/${seq}.indels.head2.vcf
GATK_INDELS_FAIL=$seq_path/${seq}/unique/${seq}.indels.FAIL.vcf

#check directory structure from previous script

## create directory structure for each sequence file
if [ ! -d "$seq_path/${seq}" ]; then
	mkdir -p $seq_path/${seq}
fi

if [ ! -d "$seq_path/${seq}/unique" ]; then
    mkdir -p $seq_path/${seq}/unique
fi  
  
#########################################################################
##   -infile(s) - SE of PE NGS reads $READ1_FILE and $READ2_file       ## 
##                                                                     ##
##                                                                     ## 
## Illumina1.8+ and Illumina old alignment with BWA for both PE and SE ##   
##             and alignment file processing with SAMtools             ##
##                                                                     ##
##      -outfiles(s) aligned sequence file $SAM and $BAM               ##
#########################################################################

if [ "$tech" == Illumina -o "$tech" == Illumina_old ]; then
    if [ "$tech" == Illumina -a "$pairing" == PE ]; then
        if [ ! -s "$SAI1" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA aln $REF_FILE $READ1_FILE -t 1 -f $SAI1 && $BWA aln $REF_FILE $READ2_FILE -t 1 -f $SAI2"
        fi
	fi
	if [ "$tech" == Illumina -a "$pairing" == SE ]; then
        if [ ! -s "$SAI1" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA aln $REF_FILE $READ1_FILE -t 1 -f $SAI1"
	    fi
	fi
	if [ "$tech" == Illumina_old -a "$pairing" == PE ]; then
        if [ ! -s "$SAI1" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA aln -I $REF_FILE $READ1_FILE -t 1 -f $SAI1 && $BWA aln -I $REF_FILE $READ2_FILE -t 1 -f $SAI2"
        fi
	fi
	if [ "$tech" == Illumina_old -a "$pairing" == SE ]; then
        if [ ! -s "$SAI1" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA aln -I $REF_FILE $READ1_FILE -t 1 -f $SAI1"
	    fi
	fi
    if [ "$pairing" == PE ]; then
	    if [ ! -s "$BAM_UNIQUE_FILE.bam" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA sampe -r '@RG\tID:${org}\tSM:${seq}\tPL:ILLUMINA\n@HD\tVN:1.0\tSO:coordinate\n@PG\tID:BWA\tPN:BWA\tVN:0.5.9' $REF_FILE $SAI1 $SAI2 $READ1_FILE $READ2_FILE > $SAM"
            log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -q 1 $SAM | $SAMTOOLS sort - $BAM_UNIQUE_FILE"
        fi
	fi
	if [ "$pairing" == SE ]; then 
	    if [ ! -s "$BAM_UNIQUE_FILE.bam" -a ! -s "$GATK_REALIGNED_BAM" ]; then
            log_eval $PBS_O_WORKDIR "$BWA samse -r '@RG\tID:${org}\tSM:${seq}\tPL:ILLUMINA\n@HD\tVN:1.0\tSO:coordinate\n@PG\tID:BWA\tPN:BWA\tVN:0.5.9' $REF_FILE $SAI1 $READ1_FILE > $SAM"
            log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -q 1 $SAM | $SAMTOOLS sort - $BAM_UNIQUE_FILE"
        fi
	fi	
    if [ ! -s "$BAM_UNIQUE_INDEX_FILE" -a ! -s "$GATK_REALIGNED_BAM" ]; then
        log_eval $PBS_O_WORKDIR "$SAMTOOLS index $BAM_UNIQUE_FILE.bam"
    fi
fi


## PGM and 454 single end

if [ "$tech" == PGM -o "$tech" == 454 -a "$pairing" == SE ]; then
    if [ ! -s "$SAM" -a ! -s "$GATK_REALIGNED_BAM" ]; then
        log_eval $PBS_O_WORKDIR "$BWA bwasw -t 1 $REF_FILE $READ1_FILE > $SAM"
	fi
	if [ ! -s "$BAM" -a ! -s "$GATK_REALIGNED_BAM" ]; then
	    log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -q 1 $SAM | $SAMTOOLS sort - $BAM"
	fi
	if [ ! -s "$BAM_UNIQUE_FILE" -a ! -s "$GATK_REALIGNED_BAM" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $ADDORREPLACEREADGROUPS SORT_ORDER=coordinate INPUT=$BAM.bam OUTPUT=$BAM_UNIQUE_FILE.bam RGID=$org RGLB=1 RGPU=1 RGPL=$tech RGSM=$seq VALIDATION_STRINGENCY=SILENT"
    fi
	if [ ! -s $BAM_INDEX_FILE -a ! -s "$GATK_REALIGNED_BAM" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $BUILDBAMINDEX INPUT=$BAM_UNIQUE_FILE.bam VALIDATION_STRINGENCY=SILENT"
    fi
	if [ -s $BAM.bam -a -s $BAM_UNIQUE_FILE.bam ]; then
		rm $BAM.bam
	fi
fi
    
	
### generated unmappped .bam file that contains the unaligned reads from the alignment using samtools
if [ ! -s "$UNMAPPED.bam" -a -s "$SAM" ]; then
    log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -f 4 -F 264 $SAM > $TMP_BAM1"
    log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -f 8 -F 260 $SAM > $TMP_BAM2"
    log_eval $PBS_O_WORKDIR "$SAMTOOLS view -h -b -S -f 12 -F 256 $SAM > $TMP_BAM3"
    log_eval $PBS_O_WORKDIR "$SAMTOOLS merge -u - $TMP_BAM1 $TMP_BAM2 $TMP_BAM3 | $SAMTOOLS sort -n - $UNMAPPED"
fi
### cleanup ###
if [ -s "$SAM" -a -s "$UNMAPPED.bam" ]; then
    rm "$SAM" "$TMP_BAM1" "$TMP_BAM2" "$TMP_BAM3"
fi

### Removal of duplicates using Picard 
if [ "$REMOVE_DUPS" == 1 ]; then
  if [ ! -s "$BAM_DUPLICATESREM_FILE" -a ! -s "$GATK_REALIGNED_BAM" ]; then ## change TMP_directory here if needed
      log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $MARKDUPLICATES I=$BAM_UNIQUE_FILE.bam O=$BAM_DUPLICATESREM_FILE REMOVE_DUPLICATES=true METRICS_FILE=$DUP_METRICS VALIDATION_STRINGENCY=LENIENT"
  fi
  if [ ! -s "$BAM_DUPLICATESREM_FILE_INDEX" -a ! -s "$GATK_REALIGNED_BAM" ]; then
      log_eval $PBS_O_WORKDIR "$SAMTOOLS index $BAM_DUPLICATESREM_FILE"
  fi
fi
if [ "$REMOVE_DUPS" == 0 ]; then
  if [ ! -s "$BAM_DUPLICATESREM_FILE" -a ! -s "$GATK_REALIGNED_BAM" ]; then ## change TMP_directory here if needed
      log_eval $PBS_O_WORKDIR "cp $BAM_UNIQUE_FILE.bam $BAM_DUPLICATESREM_FILE"
  fi
  if [ ! -s "$BAM_DUPLICATESREM_FILE_INDEX" -a ! -s "$GATK_REALIGNED_BAM" ]; then
      log_eval $PBS_O_WORKDIR "$SAMTOOLS index $BAM_DUPLICATESREM_FILE"
  fi
fi
## Realignment around poorly mapped regions ###

if [ ! -s "$GATK_TARGET_FILE" -a ! -s "$GATK_REALIGNED_BAM" ]; then
    log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T RealignerTargetCreator -R $REF_FASTA -o $GATK_TARGET_FILE -I $BAM_DUPLICATESREM_FILE"
fi
if [ ! -s "$GATK_REALIGNED_BAM" ]; then
    log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T IndelRealigner -I $BAM_DUPLICATESREM_FILE -R $REF_FASTA -targetIntervals $GATK_TARGET_FILE -o $GATK_REALIGNED_BAM"
fi
### cleanup ###
if [ -f "$BAM_DUPLICATESREM_FILE" ]; then
    rm $BAM_UNIQUE_FILE.bam $BAM_DUPLICATESREM_FILE $BAM_DUPLICATESREM_FILE_INDEX $GATK_TARGET_FILE
fi
if [ -f "$BAM_UNIQUE_INDEX_FILE" ]; then
    rm $BAM_UNIQUE_INDEX_FILE
fi
if [ -f "$BAM" ]; then
    rm $BAM
fi

### Generates a coverage analysis summary ###
if [ "$tech" == PGM ]; then
    if [ ! -s "$GATK_DEPTH_FILE" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T DepthOfCoverage -mmq 30 -R $REF_FASTA -rf BadCigar -omitBaseOutput -o $seq_path/${seq}/unique/${seq} -I $GATK_REALIGNED_BAM"
    fi
else
    if [ ! -s "$GATK_DEPTH_FILE" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T DepthOfCoverage -mmq 30 -R $REF_FASTA -omitBaseOutput -o $seq_path/${seq}/unique/${seq} -I $GATK_REALIGNED_BAM"
    fi
fi

### Call SNPS and indels ###
if [ "$tech" == PGM ]; then
    if [ ! -s "$GATK_RAW_SNPS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.snps.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -glm SNP -R $REF_FASTA -T UnifiedGenotyper -rf BadCigar -I $GATK_REALIGNED_BAM -o $GATK_RAW_SNPS_FILE -ploidy 1"
    fi
    if [ ! -s "$GATK_RAW_INDELS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.indels.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -glm INDEL -R $REF_FASTA -T UnifiedGenotyper -rf BadCigar -I $GATK_REALIGNED_BAM -o $GATK_RAW_INDELS_FILE -ploidy 1"
    fi
else
    if [ ! -s "$GATK_RAW_SNPS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.snps.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -glm SNP -R $REF_FASTA -T UnifiedGenotyper -I $GATK_REALIGNED_BAM -o $GATK_RAW_SNPS_FILE -ploidy 1"
    fi
    if [ ! -s "$GATK_RAW_INDELS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.indels.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -glm INDEL -R $REF_FASTA -T UnifiedGenotyper -I $GATK_REALIGNED_BAM -o $GATK_RAW_INDELS_FILE -ploidy 1"
    fi
fi

AVG_DEPTH=`awk /^Total/'{printf $3}' $GATK_DEPTH_FILE` 

echo -e "The Average coverage of your genome alignment calculated with GATK is $AVG_DEPTH\n"

# Process the SNPs

if [ ! -s "$GATK_FILTER_SNPS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.snps.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T VariantFiltration -R $REF_FASTA -o $GATK_FILTER_SNPS_FILE -V $GATK_RAW_SNPS_FILE --clusterSize $CLUSTER_SNP --clusterWindowSize $CLUSTER_WINDOW_SNP --filterExpression \"MLEAF < $MLEAF_SNP\" --filterName \"AFFilter\" --filterExpression \"QD < $QD_SNP\" --filterName \"QDFilter\" --filterExpression \"MQ < $MQ_SNP\" --filterName \"MQFilter\" --filterExpression \"FS > $FS_SNP\" --filterName \"FSFilter\" --filterExpression \"HaplotypeScore > $HAPLO_SNP\" --filterName \"HaplotypeScoreFilter\" --filterExpression \"QUAL < $QUAL_SNP || DP < ($AVG_DEPTH/$LOW_DEPTH) || DP > ($AVG_DEPTH*$HIGH_DEPTH)\" --filterName \"StandardFilters\" --filterExpression \"MQ0 >= 4 && '((MQ0 / (1.0 * DP))>0.1)'\" --filterName \"HARD_TO_VALIDATE\""
fi 
if [ ! -s "$GATK_PASS_SNPS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_PASS/${seq}.snps.PASS.vcf" ]; then 
        header=`grep -n "#CHROM" $GATK_FILTER_SNPS_FILE | cut -d':' -f 1`
		head -n $header $GATK_FILTER_SNPS_FILE > $GATK_SNP_HEAD
        cat $GATK_FILTER_SNPS_FILE | grep PASS | cat $GATK_SNP_HEAD - > $GATK_PASS_SNPS_FILE
fi


# Reapply filters to the raw SNPs file to generate a list of failed SNP calls

if [ ! -s "$GATK_SNP_AMB" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.snps.FAIL.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T VariantFiltration -R $REF_FASTA -o $GATK_SNP_AMB -V $GATK_RAW_SNPS_FILE --clusterSize $CLUSTER_SNP --clusterWindowSize $CLUSTER_WINDOW_SNP --filterExpression \"MLEAF < $MLEAF_SNP\" --filterName \"FAIL\" --filterExpression \"QD < $QD_SNP\" --filterName \"FAIL1\" --filterExpression \"MQ < $MQ_SNP\" --filterName \"FAIL2\" --filterExpression \"FS > $FS_SNP\" --filterName \"FAIL3\" --filterExpression \"HaplotypeScore > $HAPLO_SNP\" --filterName \"FAIL4\" --filterExpression \"QUAL < $QUAL_SNP || DP < ($AVG_DEPTH/$LOW_DEPTH) || DP > ($AVG_DEPTH*$HIGH_DEPTH)\" --filterName \"FAIL5\" --filterExpression \"MQ0 >= 4 && '((MQ0 / (1.0 * DP)) > 0.1)'\" --filterName \"FAIL6\""
fi
if [ ! -s "$GATK_SNP_FAIL" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.snps.FAIL.vcf" ]; then
            header_amb=`grep -n "#CHROM" $GATK_SNP_AMB | cut -d':' -f 1`
            head -n $header_amb $GATK_SNP_AMB > $GATK_SNP_HEAD2
            cat $GATK_SNP_AMB | grep FAIL | cat $GATK_SNP_HEAD2 - > $GATK_SNP_FAIL
fi

# process the indels

if [ ! -s "$GATK_FILTER_INDELS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.indels.PASS.vcf" ]; then
        log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T VariantFiltration -R $REF_FASTA -o $GATK_FILTER_INDELS_FILE -V $GATK_RAW_INDELS_FILE --filterExpression \"MLEAF < $MLEAF_INDEL\" --filterName \"AFFilter\" --filterExpression \"MQ < $MQ_INDEL\" --filterName \"MQFilter\" --filterExpression \"QD < $QD_INDEL\" --filterName \"QDFilter\" --filterExpression \"FS > $FS_INDEL\" --filterName \"FSFilter\" --filterExpression \"MQ0 >= 4 && '((MQ0 / (1.0 * DP)) > 0.1)'\" --filterName \"HARD_TO_VALIDATE\" --filterExpression \"QUAL < $QUAL_INDEL || DP < ($AVG_DEPTH/$LOW_DEPTH_INDEL) || DP > ($AVG_DEPTH*$HIGH_DEPTH_INDEL)\" --filterName \"QualFilter\""
fi    
if [ ! -s "$GATK_PASS_INDELS_FILE" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.indels.PASS.vcf" ]; then		
		    header_indel=`grep -n "#CHROM" $GATK_FILTER_INDELS_FILE | cut -d':' -f 1`
            head -n $header_indel $GATK_FILTER_INDELS_FILE > $GATK_INDEL_HEAD
            cat $GATK_FILTER_INDELS_FILE | grep PASS | cat $GATK_INDEL_HEAD - > $GATK_PASS_INDELS_FILE
fi
if [ ! -s "$GATK_INDELS_AMB" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.indels.FAIL.vcf" ]; then
		log_eval $PBS_O_WORKDIR "$JAVA $SET_VAR $GATK -T VariantFiltration -R $REF_FASTA -o $GATK_INDELS_AMB -V $GATK_RAW_INDELS_FILE --filterExpression \"MLEAF < $MLEAF_INDEL\" --filterName \"FAIL\" --filterExpression \"MQ < $MQ_INDEL\" --filterName \"FAIL1\" --filterExpression \"QD < $QD_INDEL\" --filterName \"FAIL2\" --filterExpression \"FS > $FS_INDEL\" --filterName \"FAIL3\" --filterExpression \"MQ0 >= 4 && '((MQ0 / (1.0 * DP)) > 0.1)'\" --filterName \"FAIL4\" --filterExpression \"QUAL < $QUAL_INDEL || DP < ($AVG_DEPTH/$LOW_DEPTH_INDEL) || DP > ($AVG_DEPTH*$HIGH_DEPTH_INDEL)\" --filterName \"FAIL5\""
fi
if [ ! -s "$GATK_INDELS_FAIL" -a ! -s "$seq_path/Outputs/SNPs_indels_FAIL/${seq}.indels.FAIL.vcf" ]; then	
		header_amb_indel=`grep -n "#CHROM" $GATK_INDELS_AMB | cut -d':' -f 1`
        head -n $header_amb_indel $GATK_INDELS_AMB > $GATK_INDELS_HEAD2
        cat $GATK_INDELS_AMB | grep FAIL | cat $GATK_INDELS_HEAD2 - > $GATK_INDELS_FAIL
fi
if [ ! -s "$BED_COV_FILE" -a ! -s "${seq_path}/BEDcov/${seq}.bedcov" ]
    then
        log_eval $PBS_O_WORKDIR "$BEDTOOLS coverage -abam $GATK_REALIGNED_BAM -b $BED_WINDOW_FILE > $BED_COV_FILE"
fi

if [ "$annotate" == yes ]; then
    if [ ! -d $seq_path/${seq}/unique/annotated ]; then
	    mkdir $seq_path/${seq}/unique/annotated
	fi
	if [ ! -s "$ANNOTATED_SNPS" ]; then
	    log_eval $PBS_O_WORKDIR/${seq}/unique/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -v $variant_genome $GATK_PASS_SNPS_FILE > $ANNOTATED_SNPS"
		mv $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_genes.txt $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_genes_SNPs.txt
		mv $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_summary.html $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_summary_SNPs.html
	fi
	if [ ! -s "$ANNOTATED_INDELS" ]; then
		log_eval $PBS_O_WORKDIR/${seq}/unique/annotated "$JAVA $SET_VAR $SNPEFF eff -no-downstream -no-intergenic -ud 100 -v $variant_genome $GATK_PASS_INDELS_FILE > $ANNOTATED_INDELS"
		mv $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_genes.txt $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_genes_indels.txt
		mv $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_summary.html $PBS_O_WORKDIR/${seq}/unique/annotated/snpEff_summary_indels.html
	fi
fi

#Moves BEDcov, SNP and indel files to output location 
if [ -s "$BED_COV_FILE" ]; then
    mv "$BED_COV_FILE" $seq_path/BEDcov
fi
if [ -s "$GATK_SNP_FAIL" ]; then
    mv "$GATK_SNP_FAIL" $seq_path/Outputs/SNPs_indels_FAIL/
fi
if [ -s "$GATK_INDELS_FAIL" ]; then 
    mv "$GATK_INDELS_FAIL" $seq_path/Outputs/SNPs_indels_FAIL/
fi
if [ -s "$GATK_PASS_SNPS_FILE" ]; then
    mv "$GATK_PASS_SNPS_FILE" $seq_path/Outputs/SNPs_indels_PASS/
fi
if [ -s "$GATK_PASS_INDELS_FILE" ]; then
    mv "$GATK_PASS_INDELS_FILE" $seq_path/Outputs/SNPs_indels_PASS/
fi
if [ -f "$SAI1" ]; then
    rm "$SAI1"
fi
if [ -f "$SAI2" ]; then
    rm "$SAI2"
fi
if [ -f "$seq_path/${seq}/unique/${seq}.sample_statistics" ]; then
    rm $seq_path/${seq}/unique/${seq}.sample_statistics $seq_path/${seq}/unique/${seq}.sample_interva* $seq_path/${seq}/unique/${seq}.sample_cumulativ*
fi

[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.snps.filter.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.snps.filter.*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.snps.AMB.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.snps.AMB.*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.snp.head.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.snp.head*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.filter.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.filter*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.indel.head.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.indel.head*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.head2.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.head*
[ -f $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.AMB.vcf ] && rm $PBS_O_WORKDIR/${seq}/unique/${seq}.indels.AMB*

sleep 20
exit 0
