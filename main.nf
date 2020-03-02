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
================================================================================
                           NF-SPANDx
                             v4.0
================================================================================

Thanks for using SPANDx!!

USAGE: nextflow run dsarov/spandx --ref <reference file>

SPANDx by default expects reads to be paired end, in the following format: STRAIN_1.fastq.gz
for the first pair and STRAIN_2.fastq.gz for the second pair.
Reads not in this format will be ignored although you can use a different read
name format by specifying the --fastq parameter

SPANDx expects at least a reference file in FASTA format.

By default all read files present in the current working directory will be processed.
Sequence files within current directory will be aligned against the reference
using BWA, SNPs and indels will be called with GATK and a SNP
matrix will be generated with GATK and VCFTools


Written by Derek Sarovich and Erin Price - University of the Sunshine Coast, Sippy Downs, Australia
Please send bug reports to dereksarovich@gmail.com
If you find SPANDx useful and use it in published work please cite - SPANDx: a
genomics pipeline for comparative analysis of large haploid whole genome
re-sequencing datasets - BMC Research Notes 2014, 7:618"

################################################################################


Input Parameter:

    --fastq      Input PE read file wildcard (default: "*_{1,2}.fastq.gz")

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

  --database     If you want to annotate the variant output then set this
                 parameter to the name of the variant file in snpEff
                 (default: false)

                 Currently, database is set to $params.database

  --phylogeny    If you would like to switch off phylogenetic reconstruction
                 and just generate a list of SNPs/indels then swith this parameter
                 to false. (default: true)

                 Currently phylogeny is set to $params.phylogeny

  --window       Default window size used in the bedcov coverage assessment
                 (default: 1kb)

                 Currently window is set to $params.window

  --assemblies   Optionally include a directory of assembled genomes in the
                 analysis. Set this parameter to 'true' if you wish to included
                 assembled genomes and place all assembled genomes in a
                 subdirectory called 'assemblies'. (default: false)

                 Currently mictures is set to $params.assemblies

  --size         ARDaP can optionally down-sample your read data to
                 run through the pipeline quicker. Set to false to skip downsampling
                 (default: 1000000)

                 Currently size is set to $params.size

  --tri_allelic  Set to true if you would like tri-allelic SNPs/indels used
                 in the phylogenetic analysis (default: false)

                 Currently tri_allelic is set to $params.tri_allelic

  --indels       Set to true if you would like indels used
                 in the phylogenetic analysis (default: true)

                 Currently indels is set to $params.indels

  --mixtures     Optionally perform within species mixtures analysis.
                 Set this parameter to 'true' if you are dealing with
                 multiple strains within the same WGS sample. (default: false)

                 Currently mixtures is set to $params.mixtures

  --structural   Set to true if you would like to identify structural variants
                 Note that this step can take a considerable amount of time if you have deep sequencing data

                 Currently structural variant assessment is set to $params.structural

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
 *  Create a bunch of indices for SPANDx
 */


ref=params.ref
snpeff_database=params.database


fastq = Channel
  .fromFilePairs("${params.fastq}", flat: true)
	.ifEmpty { exit 1, """

Input read files could not be found.
Have you included the read files in the current directory and do they have the correct naming?
With the parameters specified, SPANDx is looking for reads named ${params.fastq}.
To fix this error either rename your reads to match this formatting or specify the desired format
when initializing SPANDx e.g. --fastq "*_{1,2}_sequence.fastq.gz"

  """ }

reference_file = file(params.ref)
if( !reference_file.exists() ) {
  exit 1, """
SPANDx can't find the reference file.
It is currently looking for this file --> ${params.reference}
Please check that you have included the reference file in the current directory and rerun
"""
}

if(params.annotation) {
  //check if database has been provided
  if(params.database) {
    println "Annotation has been requested. Looking for annotation database"
    process check_and_dl_database {

    label "snpeff_dl_db"
    conda = ""
    executor 'local'

    script:
    """
    bash Check_and_DL_SnpEff_database.sh ${params.database} ${baseDir} ${ref}
    """


  }

   }  else {
     exit 1, """
     SPANDx requires a snpEff database to be specified for the annotation to work correctly
     Please use the --database flag to specific a snpEff database compatable with your
     reference genome.
     A list of available databases can be found here https://sourceforge.net/projects/snpeff/files/
     Please make sure the snpEff version matches the database version.
     """
   }


  }

//load in assemblies

if (params.assemblies) {
  assembly_loc = Channel
    .fromPath("${params.assembly_loc}", checkIfExists: true)
    .ifEmpty {"No assembled genomes will be processed"}
    .map { file ->
      def id = file.name.toString().tokenize('_').get(0)
      return tuple(id, file)
    }
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
======================================================================
      Part 1B: create synthetic reads from reference files
======================================================================
*/

if (params.assemblies) {
  process Read_synthesis {
    label "art"
    tag {"$assembly.baseName"}

    input:
    set id, file(assembly) from assemblies

    output:
    set id, file("${assembly.baseName}_1_cov.fq.gz"), file("${assembly.baseName}_2_cov.fq.gz") into (alignment_assembly)

    """
    art_illumina -i ${assembly} -p -l 150 -f 30 -m 500 -s 10 -ss HS25 -na -o ${assembly.baseName}_out
    mv ${assembly.baseName}_out1.fq ${assembly.baseName}_1_cov.fq
    mv ${assembly.baseName}_out2.fq ${assembly.baseName}_2_cov.fq
    gzip ${assembly.baseName}_1_cov.fq
    gzip ${assembly.baseName}_2_cov.fq

    """
  }
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
    if (params.size) {
            """
            seqtk sample -s 11 ${forward} $params.size | gzip - > ${id}_1_cov.fq.gz
            seqtk sample -s 11 ${reverse} $params.size | gzip - > ${id}_2_cov.fq.gz
            """
     } else {
            // Rename files if not downsampled and feed into alignment channel
            """
            mv ${forward} ${id}_1_cov.fq.gz
            mv ${reverse} ${id}_2_cov.fq.gz
            """
      }
}

/*
=======================================================================
        Part 2C: Align reads against the reference with assemblies
=======================================================================
*/
if (params.assemblies) {
  process ReferenceAlignment_assembly {

    label "spandx_alignment"
    tag {"$id"}

    input:
    file ref_index from ref_index_ch
    set id, file(forward), file(reverse) from alignment.mix(alignment_assembly) // Reads

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
} else {

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
    set id, file("${id}.bedcov")
    file("${id}.bedcov") into bedcov_files

    """
    bedtools coverage -a ${refcov} -b ${dedup_bam} > ${id}.bedcov
    """
    //mosdepth --by ${refcov} output ${dedup_bam}
    //sum_depth=\$(zcat output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    //total_chromosomes=\$(zcat output.regions.bed.gz | awk '{print \$4}' | wc -l)
    //echo "\$sum_depth/\$total_chromosomes" | bc > ${id}.depth.txt
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

  if (params.annotation) {
      process AnnotateMixture {

      label "spandx_snpeff"
      tag { "$id" }
      publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: false

      input:
      set id, file("${id}.PASS.snps.indels.mixed.vcf") from filteredMixture

      output:
      set id, file("${id}.ALL.annotated.mixture.vcf") into mixtureArdapProcessing

      //Check to see if there is a databae in the default location then run
      """
      snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} ${id}.PASS.snps.indels.mixed.vcf > ${id}.ALL.annotated.mixture.vcf
      """

      //If database isn't found then check the local directory
    //  """
    //  snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $snpeff_database ${id}.PASS.snps.indels.mixed.vcf > ${id}.ALL.annotated.mixture.vcf
    //  """
      }
  }

  if (params.structural) {
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
        if (params.annoate) {
          snpEff eff -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} \${f}.vcf > \${f}.vcf.annotated
        }
      done
      """
    }
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

  if (params.annotation) {
    process AnnotateSNPs {

      // Need to split and optimize with threads

      label "spandx_snpeff"
      tag { "$id" }
      publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: false

      input:
      set id, file(snp_pass), file(snp_fail) from filteredSNPs

      output:
      set id, file("${id}.PASS.snps.annotated.vcf") into annotatedSNPs

      //Look for the annotation in the default location
      """
      snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} $snp_pass > ${id}.PASS.snps.annotated.vcf
      """
      //if not found look in the non-default location
      //  """
      //  snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff  $snp_pass > ${id}.PASS.snps.annotated.vcf
      //  """

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

    //Look for the annotation in the default location
    """
    snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} $indel_pass > ${id}.PASS.indels.annotated.vcf
    """
    //if not found look in the non-default location
    //"""
    //snpEff eff -t -nodownload -no-downstream -no-intergenic -ud 100 -v -dataDir ${baseDir}/resources/snpeff $snpeff_database $indel_pass > ${id}.PASS.indels.annotated.vcf
    //"""
  }
 }
}

/*
===========================================================================
= This process will combine all vcf files into a master VCF file
= Clean vcf files are concatenated and converted into a matrix for phylogeny programs
=
===========================================================================
*/
process Merge_bedcov {
  label "bedcov"
  tag { "$id" }
  publishDir "./Outputs/Coverage", mode: 'copy', overwrite: false

  input:
  file("*.bedcov") from bedcov_files.collect()

  output:
  file("Bedcov_merge.txt")

  """
  bash Bedcov_merge.sh
  """
}


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
    //tag { "$id" }
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
    bash Master_vcf.sh ${reference.baseName}
    gatk VariantFiltration -R ${reference} -O out.filtered.vcf -V out.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "HaplotypeScoreFilter"
    """

  }
  if (params.annotation) {
    process snp_matrix {
      label "SNP_matrix"
      publishDir "./Outputs/Phylogeny_and_annotation", mode: 'copy', overwrite: false

      input:
      set file(filtered_vcf), file(out_vcf) from snp_matrix_ch

      output:
      file("Ortho_SNP_matrix.nex")
      file("MP_phylogeny.tre")
      file("ML_phylogeny.tre") //need to count taxa to tell this to not be expected if ntaxa is < 4
      file("All_SNPs_indels_annotated.txt")

      script:
      """
      bash SNP_matrix.sh ${snpeff_database} ${baseDir}
      """
    }
 } else {
   process snp_matrix_no_annotate {
      label "SNP_matrix"
      publishDir "./Outputs/Phylogeny", mode: 'copy', overwrite: false

      input:
      set file(filtered_vcf), file(out_vcf) from snp_matrix_ch

      output:
      file("Ortho_SNP_matrix.nex")
      file("MP_phylogeny.tre")
      file("ML_phylogeny.tre") //need to count taxa to tell this to not be expected if ntaxa is < 4


      script:
      """
      bash SNP_matrix_no_annotate.sh ${baseDir}
      """
  }
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
