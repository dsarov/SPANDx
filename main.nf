#!/usr/bin/env nextflow

/*
 *
 *  Pipeline            NF-SPANDx
 *  Version             v4.0
 *  Description         A comparative genomics pipeline
 *  Authors             Derek Sarovich, Erin Price
 *
 */

log.info """
===============================================================================
                           NF-SPANDx
                             v4.0
================================================================================

Thanks for using SPANDx!!

USAGE: nextflow run dsarov/spandx --ref <reference file>

--anotation [includes annotation parameter must match the name specified in the snpEff database] --window [Window size in base pairs for BEDcoverage module]

SPANDx by default expects reads to be paired end, in the following format: STRAIN_1.fastq.gz for the first pair and STRAIN_2.fastq.gz for the second pair.
Reads not in this format will be ignored although you can use a different read name format by specifying the --fastq parameter

SPANDx expects at least a reference file in FASTA format.

By default all read files present in the current working directory will be processed.
Sequence files within current directory will be aligned against the reference using BWA, SNPs and indels will be called with GATK and a SNP
matrix will be generated with GATK and VCFTools


Written by Derek Sarovich and Erin Price - University of the Sunshine Coast, Sippy Downs, Australia
Please send bug reports to dereksarovich@gmail.com
If you find SPANDx useful and use it in published work please cite - SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618"

#################################################################


Input Parameter:

    --fastq      Input PE read file wildcard (default: *_{1,2}.fastq.gz)

                 Currently this is set to $params.fastq

    --ref        Reference genome for alignment. Must match genome used
                 in --database (default: k96243.fasta)

                 Currently you are using $params.ref

Optional Parameters:
  --annotation   Optionally output annotated variant tables.
                 If you want to annotate the variant output then
                 set this parameter to the name of the variant file in snpEff
                 (default: false)

                 Currently annotation is set to $params.annotation

  --window       Default window size used in the bedcov coverage assessment
                 (default: 1kb)

                 Currently phylogeny is set to $params.window

  --size         ARDaP can optionally down-sample your read data to
                 run through the pipeline quicker. Set to 0 to skip downsampling
                 (default: 1000000)

                 Currently you are using $params.size

  --tri_allelic  Set to true if you would like tri-allelic SNPs/indels used
                 in the phylogenetic analysis (default: false)

                 Currently phylogeny is set to $params.tri_allelic

  --indels       Set to true if you would like indels used
                 in the phylogenetic analysis (default: false)

                 Currently phylogeny is set to $params.indels

If you want to make changes to the default `nextflow.config` file
clone the workflow into a local directory and change parameters
in `nextflow.config`:

    nextflow clone dsarov/spandx outdir/

Update to the local cache of this workflow:

    nextflow pull dsarov/spandx

==================================================================
==================================================================
"""

/*  Index Section
 *  Create a bunch of indices for ARDaP
 */


ref=params.ref
snpeff=params.annotation


fastq = Channel
  .fromFilePairs("${params.fastq}", flat: true)
	.ifEmpty { exit 1, "Input read files could not be found." }

reference_file = file(params.ref)
if( !reference_file.exists() ) {
  exit 1, "The reference file does no exist: ${params.reference}"
}

/*
======================================================================
      Part 1: create reference indices, dict files and bed files
======================================================================
*/

process IndexReference {

        label "index"

        input:
        file reference from reference_file

        output:
        file "ref.*" into ref_index_ch
        file "${reference}.fai" into ref_fai_ch1
        file "${reference.baseName}.dict" into ref_dict_ch1
        file "${reference}.bed" into refcov_ch

        """
        bwa index -a is -p ref $reference
        samtools faidx $reference
        picard CreateSequenceDictionary R=$reference O=${reference.baseName}.dict
        bedtools makewindows -g ${reference}.fai -w $params.window > ${reference}.bed
        """
}

/*
=======================================================================
Part 2: read processing, reference alignment and variant identification
=======================================================================
// Variant calling sub-workflow - basically SPANDx with a tonne of updates
// Careful here, not sure if the output overwrites the symlinks
// created by Nextflow (if input is .fq.gz) and would do weird stuff?

=======================================================================
   Part 2A: Trim reads with light quality filter and remove adapters
=======================================================================
*/
process Trimmomatic {

    label "spandx_default"
    tag {"$id"}

    input:
    set id, file(forward), file(reverse) from fastq

    output:
    set id, "${id}_1.fq.gz", "${id}_2.fq.gz" into downsample

    """
    trimmomatic PE -threads $task.cpus ${forward} ${reverse} \
    ${id}_1.fq.gz ${id}_1_u.fq.gz ${id}_2.fq.gz ${id}_2_u.fq.gz \
    ILLUMINACLIP:${baseDir}/resources/all_adapters.fa:2:30:10: \
    LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
/*
=======================================================================
              Part 2B: Downsample reads to increase speed
=======================================================================
*/
process Downsample {

    label "spandx_default"
    tag { "$id" }
    publishDir "./Clean_reads", mode: 'copy', overwrite: false

    input:
    set id, file(forward), file(reverse) from downsample

    output:
    set id, file("${id}_1_cov.fq.gz"), file("${id}_2_cov.fq.gz") into (alignment, alignmentCARD)

    script:
    if (params.size > 0) {
            """
            seqtk sample -s 11 ${forward} $params.size | gzip - > ${id}_1_cov.fq.gz
            seqtk sample -s 11 ${reverse} $params.size | gzip - > ${id}_2_cov.fq.gz
            """
     } else {
            // Rename files even if not downsampled to channel into Alignment
            """
            mv ${forward} ${id}_1_cov.fq.gz
            mv ${reverse} ${id}_2_cov.fq.gz
            """
      }
}
/*
=======================================================================
               Part 2C: Align reads against the reference
=======================================================================
*/
process ReferenceAlignment {

    label "spandx_alignment"
    tag {"$id"}

    input:
    file ref_index from ref_index_ch
    set id, file(forward), file(reverse) from alignment // Reads

    output:
    set id, file("${id}.bam"), file("${id}.bam.bai") into dup

    """
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
    -t $task.cpus ref ${forward} ${reverse} > ${id}.sam
    samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
    samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
    samtools index ${id}.bam
    """

}
/*
=======================================================================
                       Part 2D: De-duplicate bams
=======================================================================
*/
process Deduplicate {

    label "spandx_default"
    tag { "$id" }
    publishDir "./Outputs/bams", mode: 'copy', overwrite: false

    input:
    set id, file(bam_alignment), file(bam_index) from dup

    output:
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") into (averageCoverage, variantCalling, mixturePindel, variantcallingGVCF_ch)

    """
    gatk MarkDuplicates -I "${id}.bam" -O ${id}.dedup.bam --REMOVE_DUPLICATES true \
    --METRICS_FILE ${id}.dedup.txt --VALIDATION_STRINGENCY LENIENT
    samtools index ${id}.dedup.bam
    """
}
/*
=======================================================================
              Part 2E: Calculate coverage stats
=======================================================================
*/
process ReferenceCoverage {

    label "spandx_default"
    tag { "$id" }

    input:
    file refcov from refcov_ch
    set id, file(dedup_bam), file(dedup_bam_bai) from averageCoverage

    output:
    set id, file("output.per-base.bed.gz"), file("${id}.depth.txt") into (coverageData)

    """
    mosdepth --by ${refcov} output ${dedup_bam}
    sum_depth=\$(zcat output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    total_chromosomes=\$(zcat output.regions.bed.gz | awk '{print \$4}' | wc -l)
    echo "\$sum_depth/\$total_chromosomes" | bc > ${id}.depth.txt
    """
}
/*
=======================================================================
                        Part 2F: Variant identification
=======================================================================
*/
if (params.mixtures) {

  process VariantCallingMixture {

    label "spandx_gatk"
    tag { "$id" }


    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") from variantCalling

    output:
    set id, file("${id}.raw.snps.indels.mixed.vcf"), file("${id}.raw.snps.indels.mixed.vcf.idx") into mixtureFilter
    //set id, file("${id}.raw.gvcf")
	 // file("${id}.raw.gvcf") into gvcf_files
    //val true into gvcf_complete_ch

    """
    gatk HaplotypeCaller -R ${reference} --I ${id}.dedup.bam -O ${id}.raw.snps.indels.mixed.vcf
    """
  }

  process VariantFilterMixture {

    label "spandx_gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: false

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file(variants), file(variants_index) from mixtureFilter

    output:
    set id, file("${id}.PASS.snps.indels.mixed.vcf") into filteredMixture

    // Not sure if I overlooked something, but no FAIL here

    """
    gatk VariantFiltration -R ${reference} -O ${id}.snps.indels.filtered.mixed.vcf -V $variants \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "StandardFilters"

    header=`grep -n "#CHROM" ${id}.snps.indels.filtered.mixed.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.snps.indels.filtered.mixed.vcf > snp_head
		cat ${id}.snps.indels.filtered.mixed.vcf | grep PASS | cat snp_head - > ${id}.PASS.snps.indels.mixed.vcf
    """
  }

  process AnnotateMixture {

    label "spandx_snpeff"
    tag { "$id" }
    publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: false

    input:
    set id, file("${id}.PASS.snps.indels.mixed.vcf") from filteredMixture

    output:
    set id, file("${id}.ALL.annotated.mixture.vcf") into mixtureArdapProcessing

    """
    snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff ${id}.PASS.snps.indels.mixed.vcf > ${id}.ALL.annotated.mixture.vcf
    """
  }

  process PindelProcessing {

    label "spandx_pindel"
    tag { "$id" }

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    set id, file("${id}.dedup.bam"), file(alignment_index) from mixturePindel

    output:
    file("pindel.out_D.vcf") into mixtureDeletionSummary
    file("pindel.out_TD.vcf") into mixtureDuplicationSummary

    // Pindel + threads to run a bit faster
    // In the original script, there is a pindel.out_INT, here: pindel.out_INT_final

    """
    echo -e "${id}.dedup.bam\t250\tB" > pindel.bam.config
    pindel -f ${reference} -T $task.cpus -i pindel.bam.config -o pindel.out

    rm -f pindel.out_CloseEndMapped pindel.out_INT_final

    for f in pindel.out_*; do
      pindel2vcf -r ${reference} -R ${reference.baseName} -d ARDaP -p \$f -v \${f}.vcf -e 5 -is 15 -as 50000
      snpEff eff -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff \${f}.vcf > \${f}.vcf.annotated
    done
    """
  }

  process MixtureSummariesSQL {

    label "spandx_default"
    tag { "$id" }

    input:
    set id, file(variants) from mixtureArdapProcessing
    file(pindelD) from mixtureDeletionSummary
    file(pindelTD) from mixtureDuplicationSummary

    output:
    set id, file("${id}.annotated.ALL.effects") into variants_all_ch
    set id, file("${id}.Function_lost_list.txt") into function_lost_ch1, function_lost_ch2
    set id, file("${id}.deletion_summary_mix.txt") into deletion_summary_mix_ch
    set id, file("${id}.duplication_summary_mix.txt") into duplication_summary_mix_ch

    // check additional escapes in sed command

    // Use shell directive and single quotes to declare Netflow variables as !{var}
    // prevents mucking around with escaping commands in AWK, \ need to be esacaped still

    shell:

    '''
    echo 'Effects summary'

    awk '{if (match($0,"ANN=")){print substr($0,RSTART)}}' !{variants} > all.effects.tmp
    awk -F "|" '{ print $4,$10,$11,$15 }' all.effects.tmp | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//'> annotated.ALL.effects.tmp
    grep -E "#|ANN=" !{variants} > ALL.annotated.subset.vcf
    gatk VariantsToTable -V ALL.annotated.subset.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF AD -GF DP -O ALL.genotypes.subset.table
    tail -n +2 ALL.genotypes.subset.table | awk '{ print $5,$6,$7,$8 }' > ALL.genotypes.subset.table.headerless
    paste annotated.ALL.effects.tmp ALL.genotypes.subset.table.headerless > !{id}.annotated.ALL.effects

    echo 'Identification of high confidence mutations'

    grep '|HIGH|' !{variants} > ALL.func.lost
		awk '{if (match($0,"ANN=")){print substr($0,RSTART)}}' ALL.func.lost > ALL.func.lost.annotations
		awk -F "|" '{ print $4,$11,$15 }' ALL.func.lost.annotations | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//'> ALL.func.lost.annotations.tmp
		grep -E "#|\\|HIGH\\|" !{variants} > ALL.annotated.func.lost.vcf
    gatk VariantsToTable -V ALL.annotated.func.lost.vcf -F CHROM -F POS -F REF -F ALT -F TYPE -GF GT -GF AD -GF DP -O ALL.annotated.func.lost.table
		tail -n +2 ALL.annotated.func.lost.table | awk '{ print $5,$6,$7,$8 }' > ALL.annotated.func.lost.table.headerless
		paste ALL.func.lost.annotations.tmp ALL.annotated.func.lost.table.headerless > !{id}.Function_lost_list.txt

    echo 'Summary of deletions and duplications'

    grep -v '#' !{pindelD} | awk -v OFS="\t" '{ print $1,$2 }' > d.start.coords.list
		grep -v '#' !{pindelD} | gawk 'match($0, /END=([0-9]+);/,arr){ print arr[1]}' > d.end.coords.list
		grep -v '#' !{pindelD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $2 }' > mutant_depth.D
		grep -v '#' !{pindelD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $1+$2 }' > depth.D
		paste d.start.coords.list d.end.coords.list mutant_depth.D depth.D > !{id}.deletion_summary_mix.txt

    grep -v '#' !{pindelTD} | awk -v OFS="\t" '{ print $1,$2 }' > td.start.coords.list
		grep -v '#' !{pindelTD} | gawk 'match($0, /END=([0-9]+);/,arr){ print arr[1]}' > td.end.coords.list
		grep -v '#' !{pindelTD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $2 }' > mutant_depth.TD
		grep -v '#' !{pindelTD} | awk '{ print $10 }' | awk -F":" '{print $2 }' | awk -F"," '{ print $1+$2 }' > depth.TD
		paste td.start.coords.list td.end.coords.list mutant_depth.TD depth.TD > !{id}.duplication_summary_mix.txt
    '''

  }

} else {

    // Not a mixture
    //To do split GVCF calling when phylogeny isn't called

    process VariantCalling {

      label "spandx_gatk"
      tag { "$id" }
      //publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: false, pattern: '*.gvcf'

      input:
      file reference from reference_file
      file reference_fai from ref_fai_ch1
      file reference_dict from ref_dict_ch1
      set id, file(dedup_bam), file(dedup_index) from variantCalling

      output:
      set id, file("${id}.raw.snps.vcf"), file("${id}.raw.snps.vcf.idx") into snpFilter
      set id, file("${id}.raw.indels.vcf"), file("${id}.raw.indels.vcf.idx") into indelFilter
      //file("${id}.raw.gvcf") into gvcf_files
  //    val true into gvcf_complete_ch

      // v1.4 Line 261 not included yet: gatk HaplotypeCaller -R $reference -ERC GVCF --I $GATK_REALIGNED_BAM -O $GATK_RAW_VARIANTS

      """
      gatk HaplotypeCaller -R ${reference} --ploidy 1 --I ${dedup_bam} -O ${id}.raw.snps.indels.vcf
      gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.snps.vcf -select-type SNP
      gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.indels.vcf -select-type INDEL
      """
    }

  process FilterSNPs {

    label "spandx_gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: false

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file(snps), file(snps_idx) from snpFilter

    output:
    set id, file("${id}.PASS.snps.vcf"), file("${id}.FAIL.snps.vcf") into filteredSNPs

    """
    gatk VariantFiltration -R ${reference} -O ${id}.filtered.snps.vcf -V $snps \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "MLEAF < $params.MLEAF_SNP" --filter-name "AFFilter" \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "StandardFilters"

    header=`grep -n "#CHROM" ${id}.filtered.snps.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.filtered.snps.vcf > snp_head
		cat ${id}.filtered.snps.vcf | grep PASS | cat snp_head - > ${id}.PASS.snps.vcf

    gatk VariantFiltration -R ${reference} -O ${id}.failed.snps.vcf -V $snps \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "MLEAF < $params.MLEAF_SNP" --filter-name "FAIL" \
    -filter "QD < $params.QD_SNP" --filter-name "FAIL1" \
    -filter "MQ < $params.MQ_SNP" --filter-name "FAIL2" \
    -filter "FS > $params.FS_SNP" --filter-name "FAIL3" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "FAIL5"

    header=`grep -n "#CHROM" ${id}.failed.snps.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.failed.snps.vcf > snp_head
		cat ${id}.filtered.snps.vcf | grep FAIL | cat snp_head - > ${id}.FAIL.snps.vcf
    """
  }

  process FilterIndels {

    label "spandx_gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: false

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file(indels), file(indels_idx) from indelFilter

    output:
    set id, file("${id}.PASS.indels.vcf"), file("${id}.FAIL.indels.vcf") into filteredIndels

    """
    gatk VariantFiltration -R $reference -O ${id}.filtered.indels.vcf -V $indels \
    -filter "MLEAF < $params.MLEAF_INDEL" --filter-name "AFFilter" \
    -filter "QD < $params.QD_INDEL" --filter-name "QDFilter" \
    -filter "FS > $params.FS_INDEL" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_INDEL" --filter-name "QualFilter"

    header=`grep -n "#CHROM" ${id}.filtered.indels.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.filtered.indels.vcf > snp_head
		cat ${id}.filtered.indels.vcf | grep PASS | cat snp_head - > ${id}.PASS.indels.vcf

    gatk VariantFiltration -R  $reference -O ${id}.failed.indels.vcf -V $indels \
    -filter "MLEAF < $params.MLEAF_INDEL" --filter-name "FAIL" \
    -filter "MQ < $params.MQ_INDEL" --filter-name "FAIL1" \
    -filter "QD < $params.QD_INDEL" --filter-name "FAIL2" \
    -filter "FS > $params.FS_INDEL" --filter-name "FAIL3" \
    -filter "QUAL < $params.QUAL_INDEL" --filter-name "FAIL5"

    header=`grep -n "#CHROM" ${id}.failed.indels.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.failed.indels.vcf > indel_head
		cat ${id}.filtered.indels.vcf | grep FAIL | cat indel_head - > ${id}.FAIL.indels.vcf
    """
  }


  process AnnotateSNPs {

    // Need to split and optimize with threads

    label "spandx_snpeff"
    tag { "$id" }
    publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: false

    input:
    set id, file(snp_pass), file(snp_fail) from filteredSNPs

    output:
    set id, file("${id}.PASS.snps.annotated.vcf") into annotatedSNPs

    """
    snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff $snp_pass > ${id}.PASS.snps.annotated.vcf
    """
  }

  process AnnotateIndels {
    // TO DO
    // Need to split and optimize with threads

    label "spandx_snpeff"
    tag { "$id" }
    publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: false

    input:
    set id, file(indel_pass), file(indel_fail) from filteredIndels

    output:
    set id, file("${id}.PASS.indels.annotated.vcf") into annotatedIndels

    """
    snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $params.snpeff $indel_pass > ${id}.PASS.indels.annotated.vcf
    """
  }

  process VariantSummaries {

    label "spandx_default"
    tag { "$id" }

    input:
    set id, file(perbase), file(depth) from coverageData

    output:
    set id, file("${id}.deletion_summary.txt") into deletion_summary_ch
    set id, file("${id}.duplication_summary.txt") into duplication_summary_ch

    """
		echo -e "Chromosome\tStart\tEnd\tInterval" > tmp.header
		zcat $perbase | awk '\$4 ~ /^0/ { print \$1,\$2,\$3,\$3-\$2 }' > del.summary.tmp
		cat tmp.header del.summary.tmp > ${id}.deletion_summary.txt

    covdep=\$(head -n 1 $depth)
    DUP_CUTOFF=\$(echo "\$covdep*3" | bc)

    zcat $perbase | awk -v DUP_CUTOFF="\$DUP_CUTOFF" '\$4 >= DUP_CUTOFF { print \$1,\$2,\$3,\$3-\$2 }' > dup.summary.tmp

	  i=\$(head -n1 dup.summary.tmp | awk '{ print \$2 }')
	  k=\$(tail -n1 dup.summary.tmp | awk '{ print \$3 }')
	  chr=\$(head -n1 dup.summary.tmp | awk '{ print \$1 }')

	  awk -v i="\$i" -v k="\$k" -v chr="\$chr" 'BEGIN {printf "chromosome " chr " start " i " "; j=i} {if (i==\$2 || i==\$2-1 || i==\$2-2 ) {
		i=\$3;
		}
		else {
		  print "end "i " interval " i-j;
		  j=\$2;
		  i=\$3;
		  printf "chromosome " \$1 " start "j " ";
		}} END {print "end "k " interval "k-j}' < dup.summary.tmp > dup.summary.tmp1

	  sed -i 's/chromosome\\|start \\|end \\|interval //g' dup.summary.tmp1
	  echo -e "Chromosome\\tStart\\tEnd\\tInterval" > dup.summary.tmp.header
	  cat dup.summary.tmp.header dup.summary.tmp1 > ${id}.duplication_summary.txt
    """
  }

  shell:
  process VariantSummariesSQL {

    label "spandx_default"
    tag { "$id" }

    input:
    set id, file(indels) from annotatedIndels
    set id, file(snps) from annotatedSNPs

    output:
    set id, file("${id}.annotated.indel.effects") into annotated_indels_ch
    set id, file("${id}.annotated.snp.effects") into annotated_snps_ch
    set id, file("${id}.Function_lost_list.txt") into function_lost_ch1, function_lost_ch2

    shell:

    '''
    awk '{
			if (match($0,"ANN=")){print substr($0,RSTART)}
			}' !{indels} > indel.effects.tmp

		awk -F "|" '{ print $4,$10,$11,$15 }' indel.effects.tmp | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//'> !{id}.annotated.indel.effects

		awk '{
			if (match($0,"ANN=")){print substr($0,RSTART)}
			}' !{snps} > snp.effects.tmp
		awk -F "|" '{ print $4,$10,$11,$15 }' snp.effects.tmp | sed 's/c\\.//' | sed 's/p\\.//' | sed 's/n\\.//' > !{id}.annotated.snp.effects

		echo 'Identifying high consequence mutations'

		grep 'HIGH' snp.effects.tmp  | awk -F"|" '{ print $4,$11 }' >> !{id}.Function_lost_list.txt
		grep 'HIGH' indel.effects.tmp | awk -F"|" '{ print $4,$11 }' >> !{id}.Function_lost_list.txt

		sed -i 's/p\\.//' !{id}.Function_lost_list.txt
    '''

  }
}

/*
====================================================================
                              Part 3
  These processes will interrogate the SQL databases (except CARD)
  These have been split to run across different flavours of variants
                 so they can be run in parallel
=====================================================================
*/
 if (params.mixtures) {

  process SqlSnpsIndelsMix {

    label "genomic_queries"
    tag { "$id" }

    input:
    set id, file("${id}.annotated.ALL.effects") from variants_all_ch
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch1
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output_snp_indel_mix.txt") into abr_report_snp_indel_mix_ch

    script:
    """
    chmod +x ${baseDir}/bin/SQL_queries_SNP_indel_mix.sh
    SQL_queries_SNP_indel_mix.sh ${id} ${resistance_db}
    """
  }

  process SqlDeletionDuplicationMix {

    label "genomic_queries"
    tag { "$id" }

    input:
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch2
    set id, file("${id}.deletion_summary_mix.txt") from deletion_summary_mix_ch
    set id, file("${id}.duplication_summary_mix.txt") from duplication_summary_mix_ch
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output_del_dup_mix.txt") into abr_report_del_dup_mix_ch

    script:
    """
    chmod +x ${baseDir}/bin/SQL_queries_DelDupMix.sh
    SQL_queries_DelDupMix.sh ${id} ${resistance_db}
    """
  }

  process AbrReportMix {

    label "report"
    tag { "$id" }
    //publishDir "./Outputs/AbR_reports", mode: 'copy', overwrite: false

    input:
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch
    set id, file("${id}.AbR_output_del_dup_mix.txt") from abr_report_del_dup_mix_ch
    set id, file("${id}.AbR_output_snp_indel_mix.txt") from abr_report_snp_indel_mix_ch
    file("patientMetaData.csv") from patient_meta_file
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output.final.txt") into r_report_ch
    file("patientMetaData.csv") into r_report_metadata_ch
    file("patientDrugSusceptibilityData.csv") into r_report_drug_data_ch

    script:
    """
    chmod +x ${baseDir}/bin/AbR_reports.sh
    AbR_reports.sh ${id} ${resistance_db}
    """
  }
}
else {
  process SqlSnpsIndels {

    label "genomic_queries"
    tag { "$id" }

    input:
    set id, file("${id}.annotated.indel.effects") from annotated_indels_ch
    set id, file("${id}.annotated.snp.effects") from annotated_snps_ch
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch1
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output_snp_indel.txt") into abr_report_snp_indel_ch
    //Not sure if the out needs to be specific for each process or can be merged easily

    script:
    """
    chmod +x ${baseDir}/bin/SQL_queries_SNP_indel.sh
    SQL_queries_SNP_indel.sh ${id} ${resistance_db}
    """

  }

  process SqlDeletionDuplication {

    label "genomic_queries"
    tag { "$id" }

    input:
    set id, file("${id}.Function_lost_list.txt") from function_lost_ch2
    set id, file("${id}.deletion_summary.txt") from deletion_summary_ch
    set id, file("${id}.duplication_summary.txt") from duplication_summary_ch
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output_del_dup.txt") into abr_report_del_dup_ch

    script:
    """
    chmod +x ${baseDir}/bin/SQL_queries_DelDup.sh
    SQL_queries_DelDup.sh ${id} ${resistance_db}
    """
  }

  process AbrReport {

    label "report"
    tag { "$id" }
    //publishDir "./Outputs/AbR_reports", mode: 'copy', overwrite: false

    input:
    set id, file("${id}.CARD_primary_output.txt") from abr_report_card_ch
    set id, file("${id}.AbR_output_del_dup.txt") from abr_report_del_dup_ch
    set id, file("${id}.AbR_output_snp_indel.txt") from abr_report_snp_indel_ch
    file("patientMetaData.csv") from patient_meta_file
    file resistance_db from resistance_database_file

    output:
    set id, file("${id}.AbR_output.final.txt") into r_report_ch
    file("patientMetaData.csv") into r_report_metadata_ch
    file("patientDrugSusceptibilityData.csv") into r_report_drug_data_ch
    //set id, file("${id}.AbR_output.txt")

    script:
    """
    chmod +x ${baseDir}/bin/AbR_reports.sh
    AbR_reports.sh ${id} ${resistance_db}
    """
  }
}

process R_report {
  label "report"
  tag { "$id" }
  publishDir "./Outputs/AbR_reports", mode: 'copy', overwrite: true

  input:
  set id, file("${id}.AbR_output.final.txt") from r_report_ch
  file("ARDaP_logo.png") from r_report_logo_file
  file("patientMetaData.csv") from r_report_metadata_ch
  file("patientDrugSusceptibilityData.csv") from r_report_drug_data_ch
  file("sweaveTB-WGS-Micro-Report.Rnw") from sweave_report_file

  output:
  set id, file("${id}_strain.pdf")
  set id, file("${id}.AbR_output.final.txt")

  script:
  """
  chmod +x ${baseDir}/bin/Report.R
  Report.R --no-save --no-restore --args SCRIPTPATH=${baseDir} strain=${id} output_path=./
  """
}

/*
===========================================================================
= This process will combine all vcf files into a master VCF file
= Clean vcf files are concatenated and converted into a matrix for phylogeny programs
=
===========================================================================
*/

if (params.phylogeny) {

  process VariantCallingGVCF {

    label "spandx_gatk"
    tag { "$id" }
    publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: false, pattern: '*.gvcf'

    input:
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1
    set id, file("${id}.dedup.bam"), file("${id}.dedup.bam.bai") from variantcallingGVCF_ch

    output:
    //set id, file("${id}.raw.snps.indels.mixed.vcf"), file("${id}.raw.snps.indels.mixed.vcf.idx") into mixtureFilter
    set id, file("${id}.raw.gvcf")
	  file("${id}.raw.gvcf") into gvcf_files
    //val true into gvcf_complete_ch

    """
    gatk HaplotypeCaller -R ${reference} -ERC GVCF --I ${id}.dedup.bam -O ${id}.raw.gvcf
    """
  }

  process Master_vcf {
    label "master_vcf"
    tag { "id" }
    publishDir "./Outputs/Master_vcf", mode: 'copy', overwrite: false

    input:
    file("*.raw.gvcf") from gvcf_files.collect()
    file reference from reference_file
    file reference_fai from ref_fai_ch1
    file reference_dict from ref_dict_ch1

    output:
    set file("out.filtered.vcf"), file("out.vcf") into snp_matrix_ch

    script:
    """
    chmod +x ${baseDir}/bin/Master_vcf.sh
    Master_vcf.sh ${reference.baseName}
    gatk VariantFiltration -R ${reference} -O out.filtered.vcf -V out.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "HaplotypeScoreFilter"
    """

  }
  process snp_matrix {
    label "snp_matrix"
    publishDir "./Outputs/Phylogeny_and_annotation", mode: 'copy', overwrite: false

    //TO DO add additional publishDir to have annotated outputs in correct location

    input:
    set file(filtered_vcf), file(out_vcf) from snp_matrix_ch

    output:
    file("Ortho_SNP_matrix.nex")
    file("MP_phylogeny.tre")
    file("ML_phylogeny.tre") //need to count taxa to tell this to not be expected if ntaxa is < 4
    file("All_SNPs_indels_annotated.txt")

    script:
    """
    chmod +x ${baseDir}/bin/SNP_matrix.sh
    SNP_matrix.sh $params.snpeff ${baseDir}
    """
  }
}

workflow.onComplete {
	println ( workflow.success ? "\nDone! Result files are in --> ./Outputs\n \
  Antibiotic resistance reports are in --> ./Outputs/AbR_reports\n \
  If further analysis is required, bam alignments are in --> ./Outputs/bams\n \
  Phylogenetic tree and annotated merged variants are in --> ./Outputs/Phylogeny_and_annotation\n \
  Individual variant files are in --> ./Outputs/Variants/VCFs\n" \
  : "Oops .. something went wrong" )
}