#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
 *
 *  Pipeline            NF-SPANDx
 *  Version             v4.2
 *  Description         A comparative genomics pipeline
 *  Authors             Derek Sarovich, Erin Price
 *
 */

log.info """
================================================================================
                           NF-SPANDx
                             v4.3.1
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


Required Parameters:

    --ref        Reference genome for alignment. Must match genome used
                 in --database (default: k96243.fasta)

                 Currently you are using $params.ref

Optional Parameters:

  --fastq        Input PE read file wildcard (default: "*_{1,2}.fastq.gz")

                 Currently this is set to $params.fastq

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

  --single_end   Optionally use single-end Illumina data. By default SPANDx expects
                 paired-end read data. Set this variable to the directory containing
                 single end reads

                 Currently pairing is set to $params.single_end

  --window       Default window size used in the bedcov coverage assessment
                 (default: 1kb)

                 Currently window is set to $params.window

  --assemblies   Optionally include a directory of assembled genomes in the
                 analysis. Set this parameter to 'true' if you wish to included
                 assembled genomes and place all assembled genomes in a
                 subdirectory called 'assemblies'. (default: false)

                 Currently assemblies is set to $params.assemblies

  --size         SPANDx can optionally down-sample your read data to
                 run through the pipeline quicker. Set to 0 to skip downsampling
                 (default: 0). NB the number specified here refers to the number
                 of reads kept in each pair. Genome coverage will vary with genome
                 size and sequence length.

                 Currently size is set to $params.size

  --tri_allelic  Set to true if you would like tri-allelic SNPs/indels used
                 in the phylogenetic analysis (default: false).

                 Currently tri_allelic is set to $params.tri_allelic

  --indels       Set to true if you would like indels used
                 in the phylogenetic analysis (default: true).

                 Currently indels is set to $params.indels

  --mixtures     Optionally perform within species mixtures analysis.
                 Set this parameter to 'true' if you are dealing with
                 multiple strains within the same WGS sample (default: false).

                 Currently mixtures is set to $params.mixtures

  --structural   Set to true if you would like to identify structural variants
                 Note that this step can take a considerable amount of time if
                 you have deep sequencing data (default: false).

                 Currently structural variant assessment is set to $params.structural

  --notrim       Although not generally recommended to switch off, set to true
                 if you want to skip the timmomatic step (default: false).

                 Currently notrim is set to $params.notrim

  --unaligned    Optionally output unaligned reads. Useful for identifying
                 accessory genome in comparison to a reference or removing
                 unwanted contamination from raw read data (default: false).
                 Currently experimental

                 Currently unaligned is set to $params.unaligned

  --fast         Uses a different set of tools (i.e. strobealign for alignment)
                 to complete the pipeline quicker. Currently experimental

                 Currently fast is set to $params.fast

  --ont          Set to true if you have Oxford Nanopore Technologies (ONT; nanopore) data

                 Currently ont is set to $params.ont
  --ont_dir      Optionally specify a directory containing ONT data

                 Currently ont_dir set to $params.ont_dir

  --freebayes    Optionally switch to using freebayes as a variant caller instead of the GATK

                 Currently freebayes is set to $params.freebayes



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

/*
 * Improved Nextflow code to detect and load both single-end and paired-end FASTQ files
 * into separate channels.
 */
 // Initialize the paired_end_fastq channel from file pairs
 /*
 paired_end_fastq = Channel
     .fromFilePairs("${params.fastq}", size: 2, flat: true)
     .ifEmpty {
         log.info("No paired-end FASTQ files found matching pattern: ${params.fastq}")
         return Channel.empty() // Return an empty channel to avoid null issues
     }
     .map { id, files -> tuple(id, files) }
*/

     fastq = Channel
       .fromFilePairs("${params.fastq}", flat: true)
     	.ifEmpty { exit 1, """ Input read files could not be found.
     Have you included the read files in the current directory and do they have the correct naming?
     With the parameters specified, SPANDx is looking for reads named ${params.fastq}.
     To fix this error either rename your reads to match this formatting or specify the desired format
     when initializing SPANDx e.g. --fastq "*_{1,2}_sequence.fastq.gz"

     """
     }

/*
Channel
     .fromFilePairs("${params.fastq}", size: 2, flat: true)
     .ifEmpty {
         log.info("No paired-end FASTQ files found matching pattern: ${params.fastq}") }
     .map { id, files -> tuple(id, files) }
     .set { paired_end_fastq }
     log.info("Found paired end reads matching pattern: ${params.fastq}")
*/
// --- Local copies of params for channel closures ---
def run_se = params.single_end
def run_ont = params.ont
def run_asm = params.assemblies

// --- Single-End Reads ---
def logged_se = false
def single_end_files = Channel.empty() // 1. Initialize as empty
if (run_se) {
    single_end_files = Channel.fromPath("${params.single_end_dir}/*.{fastq,fq}.gz", checkIfExists: true)
}
single_end_fastq = single_end_files
    .map { file ->
        if (!logged_se) { // 4. "Log once" logic is now safely inside the map
            log.info("Found SE data. Adding to analysis")
            logged_se = true
        }
        def id = file.name.toString().tokenize('.').get(0)
        return tuple(id, file)
    }
    .ifEmpty { // 5. Log if the channel was empty (no files found or param was false)
        log.info("No single-end FASTQ files found.")
    }

// --- ONT Reads ---
def logged_ont = false
def ont_files = Channel.empty()
if (run_ont) {
    ont_files = Channel.fromPath("${params.ont_dir}/*.{fastq,fq}.gz", checkIfExists: true)
}
ont_reads = ont_files
    .map { file ->
        if (!logged_ont) {
            log.info("Found ont data. Adding to analysis")
            logged_ont = true
        }
        def id = file.name.find(/(.+?)\.(fastq|fq)\.gz$/) { it[1] }
        return tuple(id, file)
    }
    .ifEmpty {
        log.info("No ont files found or ont analysis was not requested")
    }

// --- Assembly Files ---
def logged_asm = false
def assembly_files = Channel.empty()
if (run_asm) { // This is a boolean check (true/false)
    // This now looks for uncompressed files in the 'assemblies' directory
    assembly_files = Channel.fromPath("assemblies/*.{fasta,fa,fna}", checkIfExists: true)
}
assembly_ch = assembly_files
    .map { file ->
        if (!logged_asm) {
            log.info("Found assembly data. Adding to analysis")
            logged_asm = true
        }
        // Removed .gz from the regex
        def id = file.name.find(/(.+?)\.(fasta|fa|fna)$/) { it[1] }
        return tuple(id, file)
    }
    .ifEmpty {
      log.info("No assembly files found or assembly analysis was not requested")
    }
/*
Input read files could not be found.
Please check that you have included the read files in the current directory and that they have the correct naming.
SPANDx is looking for reads named ${params.fastq}.
To fix this error, either rename your reads to match the expected formatting or specify the desired format
when initializing SPANDx e.g. --fastq "*_1_sequence.fastq.gz" for single-end or "*_{1,2}_sequence.fastq.gz" for paired-end.
*/




reference_file = file(params.ref)
reference = params.ref
reference_name = reference_file.baseName

if( !reference_file.exists() ) {
  exit 1, """
SPANDx can't find the reference file.
It is currently looking for this file --> ${params.ref}
Please check that you have included the reference file in the current directory and rerun
"""
}


process check_and_dl_database {
    label "snpeff_dl_db"

    input:
    file reference

    conda = ""
    executor 'local'

    script:
    """
    bash ${baseDir}/bin/Check_and_DL_SnpEff_database.sh ${params.database} ${baseDir} ${ref} ${params.snpeff_cache} ${params.snpeff_config}
    """
}

/*
======================================================================
      Part 1: create reference indices, dict files and bed files
======================================================================
*/

process IndexReference {


        label "index"
        conda System.getenv('CONDA_PREFIX')

        input:
        path reference_file

        output:
        path "ref.*", emit: ref_indices
        path "*.fai", emit: fai_files
        path "*.dict", emit: dict_files
        path "*.bed", emit: bed_files

        script:
        reference_path = file(params.ref)
        reference_name = reference_path.baseName
        """
        bwa index -a is -p ref $reference
        samtools faidx $reference
        picard CreateSequenceDictionary R=$reference O=${reference_name}.dict
        bedtools makewindows -g ${reference}.fai -w $params.window > ${reference}.bed
        """
}


/*
======================================================================
      Part 1B: create synthetic reads from reference files
======================================================================
*/


process Read_synthesis {

    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    tuple val(id), file(assembly_ch)

    output:
    tuple val(id), path("${id}_1_cov.fq.gz"), path("${id}_2_cov.fq.gz")

    """
    art_illumina -i ${assembly_ch} -p -l 150 -f 30 -m 500 -s 10 -ss HS25 -na -o ${id}_out
    mv ${id}_out1.fq ${id}_1_cov.fq
    mv ${id}_out2.fq ${id}_2_cov.fq
    gzip ${id}_1_cov.fq
    gzip ${id}_2_cov.fq
    """
  }


/*
=======================================================================
Part 2: read processing, reference alignment and variant identification
=======================================================================
// Variant calling sub-workflow - basically SPANDx with a tonne of updates
// Careful here, not sure if the output overwrites the symlinks
// created by Nextflow (if input is .fq.gz) and would do weird stuff?
*/


/*
=======================================================================
   Part 2A: Trim reads with light quality filter and remove adapters
=======================================================================
*/
process Trimmomatic {


    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    tuple val(id), file(forward), file(reverse)

    output:
    tuple val(id), path("${id}_1.fq.temp.gz"), path("${id}_2.fq.temp.gz")

    script:
    if (params.notrim) {
        """
        mv ${forward} ${id}_1.fq.temp.gz
        mv ${reverse} ${id}_2.fq.temp.gz
        """
    } else {
        """
        trimmomatic PE -threads 1 ${forward} ${reverse} \
        ${id}_1.fq.temp.gz ${id}_1_u.fq.gz ${id}_2.fq.temp.gz ${id}_2_u.fq.gz \
        ILLUMINACLIP:${projectDir}/resources/all_adapters.fa:2:30:10: \
        LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
        rm ${id}_1_u.fq.gz ${id}_2_u.fq.gz
        """
    }
}

// ont processing

process minimapAlign {


    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/bams", mode: 'copy', pattern: "*.dedup*", overwrite: true

    input:
    file ref_indices
    tuple val(id), file(ont_reads)
    file reference
    file bed_files

    output:
    tuple val(id), path("${id}.dedup.bam"), path("${id}.dedup.bam.bai"), emit: dedup_bams
    tuple val(id), path("${id}.depth.txt"), emit:depth_files

    script:
    """
    minimap2 -ax map-ont -t ${task.cpus} ${reference} ${ont_reads} | samtools sort -o ${id}.dedup.bam
    samtools index ${id}.dedup.bam
    mosdepth --by ${ref}.bed ${id}.output ${id}.dedup.bam
    sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
    echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
    """
}

/*
=======================================================================
              Part 2B: Downsample reads to increase speed
=======================================================================
*/
  process Downsample {


      label "spandx_default"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }
    // publishDir "./Clean_reads", mode: 'copy', overwrite: false

      input:
      tuple val(id), file(forward), file(reverse)

      output:
      tuple val(id), path("${id}_1_cov.fq.gz"), path("${id}_2_cov.fq.gz")

      script:
      if (params.size > 0) {
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
               Part 2C: Align reads against the reference
=======================================================================
*/

process ReferenceAlignment_assembly {


  label "spandx_alignment"
  conda System.getenv('CONDA_PREFIX')
  tag { "$id" }
  publishDir "./Outputs/bams", mode: 'copy', pattern: "*.dedup*", overwrite: true

  input:
  file ref_indices
  tuple val(id), file(forward), file(reverse)
  file reference
  file bed_files

  output:
  tuple val(id), path("${id}.dedup.bam"), path("${id}.dedup.bam.bai"), emit: dedup_bams
  tuple val(id), path("${id}.depth.txt"), emit:depth_files

  script:
  if(params.fast) {
  """
  strobealign --rg-id ID:${id} --rg SM:${id} --rg PL:ILLUMINA ${ref} ${forward} ${reverse} > ${id}.sam
  samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
  samtools sort -@ 1 -o ${id}.dedup.bam ${id}.bam_tmp
  samtools index ${id}.dedup.bam
  rm ${id}.sam ${id}.bam_tmp
  mosdepth --by ${ref}.bed ${id}.output ${id}.dedup.bam
  sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
  total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
  echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
  """
  } else {
  """
  bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
  -t $task.cpus ref ${forward} ${reverse} > ${id}.sam
  samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
  samtools sort -@ 1 -o ${id}.dedup.bam ${id}.bam_tmp
  samtools index ${id}.dedup.bam
  rm ${id}.sam ${id}.bam_tmp
  mosdepth --by ${ref}.bed ${id}.output ${id}.dedup.bam
  sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
  total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
  echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
  """
  }
}

process ReferenceAlignment {


    label "spandx_alignment"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    file ref_indices
    tuple val(id), file(forward), file(reverse)
    file reference
    file bed_files

    output:
    tuple val(id), path("${id}.bam"), path("${id}.bam.bai"), emit:bams
    tuple val(id), path("${id}.depth.txt"), emit:depth_files

    script:
    if(params.fast) {
    """
    strobealign  --rg-id ID:${id} --rg SM:${id} --rg PL:ILLUMINA ${ref} ${forward} ${reverse} > ${id}.sam
    samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
    samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
    samtools index ${id}.bam
    rm ${id}.sam ${id}.bam_tmp
    mosdepth --by ${ref}.bed ${id}.output ${id}.bam
    sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
    echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
    """
    } else {
    """
    bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
    -t $task.cpus ref ${forward} ${reverse} > ${id}.sam
    samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
    samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
    samtools index ${id}.bam
    rm ${id}.sam ${id}.bam_tmp
    mosdepth --by ${ref}.bed ${id}.output ${id}.bam
    sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
    total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
    echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
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
process Trimmomatic_SE {

    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    tuple val(id), path(forward)

    output:
    tuple val(id), path("${id}_1.fq.gz")

    script:
    if (params.notrim) {
        """
        mv $forward ${id}_1.fq.gz
        """
    } else {
        """
        trimmomatic SE -threads 1 $forward ${id}_1.fq.gz \
        ILLUMINACLIP:${baseDir}/resources/all_adapters.fa:2:30:10: \
        LEADING:10 TRAILING:10 SLIDINGWINDOW:4:15 MINLEN:36
        """
    }
}
  /*
  =======================================================================
                Part 2B: Downsample reads to increase speed
  =======================================================================
  */
process Downsample_SE {

    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    tuple val(id), path(forward)

    output:
    tuple val(id), path("${id}_1_cov.fq.gz")

    script:
    if (params.size > 0) {
        """
        seqtk sample -s 11 $forward $params.size | gzip - > ${id}_1_cov.fq.gz
        """
    } else {
        """
        mv $forward ${id}_1_cov.fq.gz
        """
    }
}
  /*
  =======================================================================
                 Part 2C: Align reads against the reference
  =======================================================================
  */
process SE_reference_alignment {


      label "spandx_alignment"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }

      input:
      path ref_indices
      tuple val(id), path(forward)
      path reference
      path bed_files

      output:
      tuple val(id), path("${id}.bam"), path("${id}.bam.bai")
      tuple val(id), path("${id}.depth.txt")

	  script:
      """
      bwa mem -R '@RG\\tID:${params.org}\\tSM:${id}\\tPL:ILLUMINA' -a \
      -t $task.cpus ref ${forward} > ${id}.sam
      samtools view -h -b -@ 1 -q 1 -o ${id}.bam_tmp ${id}.sam
      samtools sort -@ 1 -o ${id}.bam ${id}.bam_tmp
      samtools index ${id}.bam
      rm ${id}.sam ${id}.bam_tmp
      mosdepth --by ${ref}.bed ${id}.output ${id}.bam
      sum_depth=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | awk '{s+=\$1}END{print s}')
      total_chromosomes=\$(zcat ${id}.output.regions.bed.gz | awk '{print \$4}' | wc -l)
      echo "\$sum_depth \$total_chromosomes" | awk '{print int(\$1/\$2)}' > "${id}.depth.txt"
      """
}

/*
=======================================================================
                       Part 2D: De-duplicate bams
=======================================================================
*/
process Deduplicate {


    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/bams", mode: 'copy', pattern: "*.dedup*", overwrite: true

    if (params.unaligned) {
    publishDir "./Outputs/bams", mode: 'copy', pattern: "*.unmapped.bam", overwrite: true
    publishDir "./Outputs/unmapped_reads", mode: 'copy', pattern: "*fastq.gz", overwrite: true
    }

    input:
    tuple val(id), path(bam_file), path(bai_file)

    output:
    tuple val(id), path("${id}.dedup.bam"), path("${id}.dedup.bam.bai"), emit: dedup_bams

    if (params.unaligned) {
    val(id), file("${id}_unmapped_1_sequence.fastq.gz"), file("${id}_unmapped_2_sequence.fastq.gz")
    }

    script:
    if (params.unaligned) {
    """
    gatk --java-options -Xmx${task.memory.toString().replaceAll(/[\sB]/,'')} MarkDuplicates -I "${id}.bam" -O ${id}.dedup.bam --REMOVE_DUPLICATES true \
    --METRICS_FILE ${id}.dedup.txt --VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER coordinate
    samtools index ${id}.dedup.bam
    samtools view -b -F 2 ${id}.dedup.bam | samtools sort -n - -o ${id}.unmapped.bam
    bamToFastq -i ${id}.unmapped.bam -fq ${id}_unmapped_1_sequence.fastq -fq2 ${id}_unmapped_2_sequence.fastq
    gzip ${id}_unmapped_1_sequence.fastq
    gzip ${id}_unmapped_2_sequence.fastq
    """
    } else {
    """
    gatk --java-options -Xmx${task.memory.toString().replaceAll(/[\sB]/,'')} MarkDuplicates -I "${id}.bam" -O ${id}.dedup.bam --REMOVE_DUPLICATES true \
    --METRICS_FILE ${id}.dedup.txt --VALIDATION_STRINGENCY LENIENT --ASSUME_SORT_ORDER coordinate
    samtools index ${id}.dedup.bam
    """
  }
}
/*
=======================================================================
              Part 2E: Calculate coverage stats
=======================================================================
*/
process ReferenceCoverage {


    label "spandx_default"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }

    input:
    file ref_indices
    tuple val(id), path(dedup_bam), path(dedup_bam_bai)

    output:
    tuple val(id), path("*.bedcov"), emit: bedcov_files

    script:
    """
    bedtools coverage -sorted -a ${ref_indices[0]} -b ${dedup_bam} > ${id}.bedcov
    """
}

process Merge_bedcov {


  label "bedcov"
  conda System.getenv('CONDA_PREFIX')
// tag { "$id" }
  publishDir "./Outputs/Coverage", mode: 'copy', overwrite: true

  input:
  file bedcov_files

  output:
  path "Bedcov_merge.txt"

  """
  bash ${baseDir}/bin/Bedcov_merge.sh
  """
}
/*
=======================================================================
                        Part 2F: Variant identification
=======================================================================
*/
process VariantCallingMixture {


    label "spandx_gatk"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: true, pattern: '*.gvcf'


    input:
    file reference
    file fai_files
    file dict_files
    tuple val(id), path(bam_file), path(bai_file)

    output:
    tuple val(id), path("${id}.raw.snps.indels.mixed.gvcf"), emit: gvcf_files
    tuple val(id), path("${id}.raw.snps.indels.mixed.vcf"), path("${id}.raw.snps.indels.mixed.vcf.idx"), emit: vcf_files

    script:
    if (params.freebayes) {
    """
    freebayes --gvcf -f ${reference} ${id}.dedup.bam > ${id}.raw.snps.indels.mixed.gvcf
    bcftools view ${id}.raw.snps.indels.mixed.gvcf | awk '{gsub(/<\\*>/, "<NON_REF>"); print}' > converted_output.g.vcf
    mv converted_output.g.vcf ${id}.raw.snps.indels.mixed.gvcf
    gatk SelectVariants -V ${id}.raw.snps.indels.mixed.gvcf -O ${id}.raw.snps.indels.mixed.vcf --exclude-non-variants
    """
    } else {
    """
     gatk HaplotypeCaller -ERC GVCF -R ${reference} --I ${id}.dedup.bam -O ${id}.raw.snps.indels.mixed.gvcf
     gatk SelectVariants -V ${id}.raw.snps.indels.mixed.gvcf -O ${id}.raw.snps.indels.mixed.vcf --exclude-non-variants
    """
    }
}

process VariantCallingMixture_Clair3 {

    conda "$baseDir/env_clair3.yaml"
    label "spandx_clair3"
    tag { "$id" }
    publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: true, pattern: '*.gvcf'
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: true, pattern: '*.vcf'

    input:
    file reference
    file fai_files
    file dict_files
    tuple val(id), path(bam_file), path(bai_file)

    output:
    tuple val(id), path("${id}.raw.snps.indels.mixed.gvcf"), emit: gvcf_files
    tuple val(id), path("${id}.raw.snps.indels.mixed.vcf"), path("${id}.raw.snps.indels.mixed.vcf.idx"), emit: vcf_files
    tuple val(id), path("${id}.raw.snps.vcf"), path("${id}.raw.snps.vcf.idx"), emit: snps_output
    tuple val(id), path("${id}.raw.indels.vcf"), path("${id}.raw.indels.vcf.idx"), emit: indels_output

    script:
    """
	CWD=\$(pwd)
    run_clair3.sh --bam_fn ${id}.dedup.bam --ref_fn \${CWD}/${reference} --threads $task.cpus \
    --model_path "${baseDir}"/resources/clair3_models/r941_prom_sup_g5014 --output ./ -p ont --fast_mode \
    --gvcf --enable_long_indel --sample_name="${id}" --include_all_ctgs
    gunzip *.gz
    mv merge_output.vcf "${id}.raw.snps.indels.mixed.vcf"
    mv merge_output.gvcf "${id}.raw.snps.indels.mixed.gvcf"
    gatk IndexFeatureFile -I "${id}.raw.snps.indels.mixed.vcf"
    gatk IndexFeatureFile -I "${id}.raw.snps.indels.mixed.gvcf"

    gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.mixed.vcf -O ${id}.raw.snps.vcf -select-type SNP
    gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.mixed.vcf -O ${id}.raw.indels.vcf -select-type INDEL
    """
}

process VariantCalling_Clair3 {

    conda "$baseDir/env_clair3.yaml"
    label "spandx_clair3"
    tag { "$id" }
    publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: true, pattern: '*.gvcf'
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: true, pattern: '*.vcf'

    input:
    file reference
    file fai_files
    file dict_files
    tuple val(id), path(bam_file), path(bai_file)

    output:
    tuple val(id), path("${id}.raw.snps.indels.mixed.gvcf"), emit: gvcf_files
    tuple val(id), path("${id}.raw.snps.indels.mixed.vcf"), path("${id}.raw.snps.indels.mixed.vcf.idx"), emit: vcf_files
    tuple val(id), path("${id}.raw.snps.vcf"), path("${id}.raw.snps.vcf.idx"), emit: snps_output
    tuple val(id), path("${id}.raw.indels.vcf"), path("${id}.raw.indels.vcf.idx"), emit: indels_output

    script:
    """
    run_clair3.sh --bam_fn ${task.workDir}/${id}.dedup.bam --ref_fn ${task.workDir}/${reference} --threads $task.cpus \
    --model_path "${baseDir}"/resources/clair3_models/r941_prom_sup_g5014 --output ./ -p ont --fast_mode \
    --gvcf --enable_long_indel --sample_name="${id}" --include_all_ctgs
    gunzip *.gz
    mv merge_output.vcf "${id}.raw.snps.indels.mixed.vcf"
    mv merge_output.gvcf "${id}.raw.snps.indels.mixed.gvcf"
    gatk IndexFeatureFile -I "${id}.raw.snps.indels.mixed.vcf"
    gatk IndexFeatureFile -I "${id}.raw.snps.indels.mixed.gvcf"

    gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.mixed.vcf -O ${id}.raw.snps.vcf -select-type SNP
    gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.mixed.vcf -O ${id}.raw.indels.vcf -select-type INDEL
    """
}

process VariantFilterMixture {


    label "spandx_gatk"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: true

    input:
    file reference
    file reference_fai
    file reference_dict
    tuple val(id), file(variants), file(variants_index)

    output:
    tuple val(id), path("${id}.PASS.snps.indels.mixed.vcf"), path("${id}.FAIL.snps.indels.mixed.vcf")

    script:
    """
    gatk VariantFiltration -R ${reference} -O ${id}.snps.indels.filtered.mixed.vcf -V $variants \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "StandardFilters"

    header=`grep -n "#CHROM" ${id}.snps.indels.filtered.mixed.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.snps.indels.filtered.mixed.vcf > snp_head
		cat ${id}.snps.indels.filtered.mixed.vcf | awk -F'\t' '\$7 == "PASS" {print}' | cat snp_head - > ${id}.PASS.snps.indels.mixed.vcf
    cat ${id}.snps.indels.filtered.mixed.vcf | awk -F'\t' '\$7 != "PASS" {print}' | cat snp_head - > ${id}.FAIL.snps.indels.mixed.vcf
    """
}

process AnnotateMixture {


      label "spandx_snpeff"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }
      publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: true

      input:
      tuple val(id), path("${id}.PASS.snps.indels.mixed.vcf"), path("${id}.FAIL.snps.indels.mixed.vcf")

      output:
      tuple val(id), path("${id}.ALL.annotated.mixture.vcf")

      //Check to see if there is a databae in the default location then run
      script:
      """
      snpEff eff -dataDir ${params.snpeff_cache} -c ${params.snpeff_config} -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} ${id}.PASS.snps.indels.mixed.vcf > ${id}.ALL.annotated.mixture.vcf
      """
}



// TO DO - needs to be updated with Delly
process PindelProcessing {


      label "spandx_pindel"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }

      input:
      path reference
      path reference_fai
      tuple val(id), path("${id}.dedup.bam"), path(alignment_index)

      output:
      path("pindel.out_D.vcf")
      path("pindel.out_TD.vcf")


      script:
      """
      echo -e "${id}.dedup.bam\t250\tB" > pindel.bam.config
      pindel -f ${reference} -T $task.cpus -i pindel.bam.config -o pindel.out

      rm -f pindel.out_CloseEndMapped pindel.out_INT_final

      for f in pindel.out_*; do
        pindel2vcf -r ${reference} -R ${reference.baseName} -d ARDaP -p \$f -v \${f}.vcf -e 5 -is 15 -as 50000
      if [ "${params.annotate}" == "true" ]; then
            snpEff eff -dataDir ${params.snpeff_cache} -no-downstream -no-intergenic -ud 100 -v ${params.snpeff_database} \${f}.vcf > \${f}.vcf.annotated
      fi
      done
      """
}

    //Not a mixture
process VariantCalling {


      label "spandx_gatk"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }
      publishDir "./Outputs/Variants/GVCFs", mode: 'copy', overwrite: false, pattern: '*.gvcf'

      input:
      file reference
      file fai_files
      file dict_files
      tuple val(id), path(dedup_bam), path(bai_file)

      output:
      tuple val(id), path("${id}.raw.snps.vcf"), path("${id}.raw.snps.vcf.idx"), emit: snps_output
      tuple val(id), path("${id}.raw.indels.vcf"), path("${id}.raw.indels.vcf.idx"), emit: indels_output
      tuple val(id), path("${id}.raw.snps.indels.gvcf"), emit: gvcf_files

      script:
      if (params.freebayes) {
      """
      freebayes --gvcf -f ${reference} ${id}.dedup.bam > ${id}.raw.snps.indels.gvcf
      bcftools view ${id}.raw.snps.indels.mixed.gvcf | awk '{gsub(/<\\*>/, "<NON_REF>"); print}' > converted_output.g.vcf
      mv converted_output.g.vcf ${id}.raw.snps.indels.mixed.gvcf
      gatk SelectVariants -V ${id}.raw.snps.indels.gvcf -O ${id}.raw.snps.indels.vcf --exclude-non-variants
      gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.snps.vcf -select-type SNP
      gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.vcf -O ${id}.raw.indels.vcf -select-type INDEL
      """
      } else {
      """
       gatk HaplotypeCaller -R ${reference} -ERC GVCF --ploidy 1 --I ${dedup_bam} -O ${id}.raw.snps.indels.gvcf
       ## gatk SelectVariants -V ${id}.raw.snps.indels.gvcf -O ${id}.raw.snps.indels.vcf --exclude-non-variants
       gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.gvcf -O ${id}.raw.snps.vcf -select-type SNP --ignore-non-ref-in-types
       gatk SelectVariants -R ${reference} -V ${id}.raw.snps.indels.gvcf -O ${id}.raw.indels.vcf -select-type INDEL --ignore-non-ref-in-types
      """
      }
}

process FilterSNPs {


    label "spandx_gatk"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: true

    input:
    file reference
    file reference_fai
    file reference_dict
    tuple val(id), path(snps_vcf), path(snp_vcf_idx)

    output:
    tuple val(id), path("${id}.PASS.snps.vcf"), path("${id}.FAIL.snps.vcf")

    script:
    """
    gatk VariantFiltration -R ${reference} -O ${id}.filtered.snps.vcf -V ${id}.raw.snps.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "MLEAF < $params.MLEAF_SNP" --filter-name "AFFilter" \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "StandardFilters"

    header=`grep -n "#CHROM" ${id}.filtered.snps.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.filtered.snps.vcf > snp_head
		cat ${id}.filtered.snps.vcf | grep PASS | cat snp_head - > ${id}.PASS.snps.vcf

    gatk VariantFiltration -R ${reference} -O ${id}.failed.snps.vcf -V ${id}.raw.snps.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "MLEAF < $params.MLEAF_SNP" --filter-name "FAIL" \
    -filter "QD < $params.QD_SNP" --filter-name "FAIL1" \
    -filter "MQ < $params.MQ_SNP" --filter-name "FAIL2" \
    -filter "FS > $params.FS_SNP" --filter-name "FAIL3" \
    -filter "QUAL < $params.QUAL_SNP" --filter-name "FAIL5"

    header=`grep -n "#CHROM" ${id}.failed.snps.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.failed.snps.vcf > snp_head
		cat ${id}.failed.snps.vcf | grep FAIL | cat snp_head - > ${id}.FAIL.snps.vcf
    """
}

process FilterIndels {


    label "spandx_gatk"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/Variants/VCFs", mode: 'copy', overwrite: true

    input:
    file reference
    file reference_fai
    file reference_dict
    tuple val(id), path(indel_vcf), path(indel_vcf_idx)

    output:
    tuple val(id), path("${id}.PASS.indels.vcf"), path("${id}.FAIL.indels.vcf")

    script:
    """
    gatk VariantFiltration -R $reference -O ${id}.filtered.indels.vcf -V ${id}.raw.indels.vcf \
    -filter "MLEAF < $params.MLEAF_INDEL" --filter-name "AFFilter" \
    -filter "QD < $params.QD_INDEL" --filter-name "QDFilter" \
    -filter "FS > $params.FS_INDEL" --filter-name "FSFilter" \
    -filter "QUAL < $params.QUAL_INDEL" --filter-name "QualFilter"

    header=`grep -n "#CHROM" ${id}.filtered.indels.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.filtered.indels.vcf > snp_head
		cat ${id}.filtered.indels.vcf | grep PASS | cat snp_head - > ${id}.PASS.indels.vcf

    gatk VariantFiltration -R  $reference -O ${id}.failed.indels.vcf -V ${id}.raw.indels.vcf \
    -filter "MLEAF < $params.MLEAF_INDEL" --filter-name "FAIL" \
    -filter "MQ < $params.MQ_INDEL" --filter-name "FAIL1" \
    -filter "QD < $params.QD_INDEL" --filter-name "FAIL2" \
    -filter "FS > $params.FS_INDEL" --filter-name "FAIL3" \
    -filter "QUAL < $params.QUAL_INDEL" --filter-name "FAIL5"

    header=`grep -n "#CHROM" ${id}.failed.indels.vcf | cut -d':' -f 1`
		head -n "\$header" ${id}.failed.indels.vcf > indel_head
		cat ${id}.failed.indels.vcf | grep FAIL | cat indel_head - > ${id}.FAIL.indels.vcf
    """
}

process AnnotateSNPs {


      label "spandx_snpeff"
      conda System.getenv('CONDA_PREFIX')
      tag { "$id" }
      publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: true

      input:
      tuple val(id), path(snp_pass), path(snp_fail)

      output:
      tuple val(id), path("${id}.PASS.snps.annotated.vcf")

      """
      snpEff eff -dataDir ${params.snpeff_cache} -c ${params.snpeff_config} -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} $snp_pass > ${id}.PASS.snps.annotated.vcf
      """

}


process AnnotateIndels {


    label "spandx_snpeff"
    conda System.getenv('CONDA_PREFIX')
    tag { "$id" }
    publishDir "./Outputs/Variants/Annotated", mode: 'copy', overwrite: true

    input:
    tuple val(id), path(indel_pass), path(indel_fail)

    output:
    tuple val(id), path("${id}.PASS.indels.annotated.vcf")

    //Look for the annotation in the default location
    """
    snpEff eff -dataDir ${params.snpeff_cache} -c ${params.snpeff_config} -nodownload -no-downstream -no-intergenic -ud 100 -v ${snpeff_database} $indel_pass > ${id}.PASS.indels.annotated.vcf
    """

}


/*
===========================================================================
= This process will combine all vcf files into a master VCF file
= Clean vcf files are concatenated and converted into a matrix for phylogeny programs
=
===========================================================================
*/

process Master_vcf {


    label "master_vcf"
    conda System.getenv('CONDA_PREFIX')
    publishDir "./Outputs/Master_vcf", mode: 'copy', overwrite: true

    input:
    file gvcfs
    file reference_file
    file reference_fai
    file reference_dict

    output:
    tuple path("out.filtered.vcf"), path("out.vcf")

    script:
    """
    bash ${baseDir}/bin/Master_vcf.sh ${reference_name}
    gatk VariantFiltration -R ${reference} -O out.filtered.vcf -V out.vcf \
    --cluster-size $params.CLUSTER_SNP -window $params.CLUSTER_WINDOW_SNP \
    -filter "QD < $params.QD_SNP" --filter-name "QDFilter" \
    -filter "MQ < $params.MQ_SNP" --filter-name "MQFilter" \
    -filter "FS > $params.FS_SNP" --filter-name "HaplotypeScoreFilter"
    """

}

process snp_matrix {

      label "SNP_matrix"
      conda System.getenv('CONDA_PREFIX')
      publishDir "./Outputs/Phylogeny_and_annotation", mode: 'copy', overwrite: true
      publishDir "./Outputs/Phylogeny_and_annotation", mode: 'copy', overwrite: true, pattern: '*.nex'

      input:
      path(filtered_vcf)
      file depthFiles

      output:
      path("Ortho_SNP_matrix.nex")
      path("MP_phylogeny.tre")
      path("ML_phylogeny.tre")
      path("All_SNPs_indels_annotated.txt")
	  path("All_SNPs_indels_annotated_mixtures.txt")
      path("indel_matrix.nex")
      path("indel_SNP_matrix.nex")
      path("QC_metrics_summary.tsv")
      path("indel_differences_matrix.tsv")
      path("snp_differences_matrix.tsv")
      path("merged_snp_indel_matrix.tsv")

      script:
      def cmd = ""
      if (params.mixtures) {
          cmd = """
          bash ${baseDir}/bin/SNP_matrix.sh ${snpeff_database} ${baseDir} ${params.snpeff_cache} ${params.snpeff_config}
		  python ${baseDir}/bin/process_vcf_mixtures.py > "All_SNPs_indels_annotated_mixtures.txt"
          bash ${baseDir}/bin/Summary.sh ${ref} ${baseDir}
          """
      } else {
          cmd = """
          bash ${baseDir}/bin/SNP_matrix.sh ${snpeff_database} ${baseDir} ${params.snpeff_cache} ${params.snpeff_config}
		  python ${baseDir}/bin/process_vcf_mixtures.py > "All_SNPs_indels_annotated_mixtures.txt"
          bash ${baseDir}/bin/Summary_no_mixtures.sh ${ref} ${baseDir}
          """
      }

      """
      $cmd
      """

}

process snp_matrix_no_annotate {

      label "SNP_matrix"
      conda System.getenv('CONDA_PREFIX')
      publishDir "./Outputs/Phylogeny", mode: 'copy', overwrite: true

      input:
      path(filtered_vcf) //, path(out_vcf)

      output:
      path("Ortho_SNP_matrix.nex")
      path("MP_phylogeny.tre")
      path("ML_phylogeny.tre")
	    path("QC_metrics_summary.tsv")
      path("snp_differences_matrix.tsv")



      script:
      """
      bash ${baseDir}/bin/SNP_matrix_no_annotate.sh ${baseDir}
	    bash ${baseDir}/bin/Summary_no_mixtures.sh ${ref} ${baseDir}
      """
}


workflow {

    // --- 1. SETUP AND INDEXING ---
    if (params.annotation) {
        if (params.database) {
            println "Annotation has been requested. Looking for annotation database"
            check_and_dl_database(reference_file)
        } else {
            exit 1, """
            SPANDx requires a snpEff database to be specified for the annotation to work correctly.
            Please use the --database flag to specify a snpEff database compatible with your
            reference genome.
            """
        }
    }
    indexResults = IndexReference(reference_file)


    // --- 2. PROCESS DATA FROM EACH SOURCE ---

    // Paired-End Illumina Data (Baseline)
    pe_trimmed = Trimmomatic(fastq)
    pe_downsampled = Downsample(pe_trimmed)
    pe_alignment_out = ReferenceAlignment(indexResults.ref_indices, pe_downsampled, reference_file, indexResults.bed_files)
    pe_dedup_out = Deduplicate(pe_alignment_out.bams)

    // Assembled Genomes Data (Optional)
    def assembly_out
    if (params.assemblies) {
        synthetic_reads = Read_synthesis(assembly_ch)
        assembly_out = ReferenceAlignment_assembly(indexResults.ref_indices, synthetic_reads, reference_file, indexResults.bed_files)
    } else {
        // If no assemblies, create a placeholder map with empty channels
        assembly_out = [ dedup_bams: Channel.empty(), depth_files: Channel.empty() ]
    }

    // ONT Data (Optional)
    def ont_out
    if (params.ont) {
        ont_out = minimapAlign(indexResults.ref_indices, ont_reads, reference_file, indexResults.bed_files)
    } else {
        // If no ONT data, create a placeholder map with empty channels
        ont_out = [ dedup_bams: Channel.empty(), depth_files: Channel.empty() ]
    }

    // --- 3. COMBINE RESULTS FROM ALL SOURCES ---
    all_deduplicated_bams = pe_dedup_out.dedup_bams.mix(assembly_out.dedup_bams, ont_out.dedup_bams)
    all_depth_files = pe_alignment_out.depth_files.mix(assembly_out.depth_files, ont_out.depth_files)


    // --- 4. DOWNSTREAM ANALYSIS ---
    bedcov_files = ReferenceCoverage(indexResults.bed_files, all_deduplicated_bams).collect()
    merged_bedcov = Merge_bedcov(bedcov_files)

    // Variant Calling and Filtering for Illumina/Assembly data
    def illumina_variants
    if (params.mixtures) {
        // Combine PE and Assembly bams for mixture calling
        illumina_bams = pe_dedup_out.dedup_bams.mix(assembly_out.dedup_bams)
        illumina_variants = VariantCallingMixture(reference_file, indexResults.fai_files, indexResults.dict_files, illumina_bams)
    } else {
        illumina_bams = pe_dedup_out.dedup_bams.mix(assembly_out.dedup_bams)
        illumina_variants = VariantCalling(reference_file, indexResults.fai_files, indexResults.dict_files, illumina_bams)
    }

    // Variant Calling for ONT data
    def ont_variants
    if (params.ont) {
        if (params.mixtures) {
            ont_variants = VariantCallingMixture_Clair3(reference_file, indexResults.fai_files, indexResults.dict_files, ont_out.dedup_bams)
        } else {
            ont_variants = VariantCalling_Clair3(reference_file, indexResults.fai_files, indexResults.dict_files, ont_out.dedup_bams)
        }
    } else {
        // Placeholder for ONT variants if not used
        ont_variants = [ gvcf_files: Channel.empty(), vcf_files: Channel.empty(), snps_output: Channel.empty(), indels_output: Channel.empty() ]
    }

    // Filter and Annotate based on mixture vs. standard mode
    if (params.mixtures) {
        all_vcfs = illumina_variants.vcf_files.mix(ont_variants.vcf_files)
        filtered_mixed_variants = VariantFilterMixture(reference_file, indexResults.fai_files, indexResults.dict_files, all_vcfs)
        if (params.annotation) {
            annotated_mixed_variants = AnnotateMixture(filtered_mixed_variants)
        }
    } else {
                // --- THIS IS THE FIX ---
                // Now we mix the split outputs from both Illumina and ONT
                all_snps = illumina_variants.snps_output.mix(ont_variants.snps_output)
                all_indels = illumina_variants.indels_output.mix(ont_variants.indels_output)

                // The rest of the pipeline now receives all variants
                filtered_snps = FilterSNPs(reference_file, indexResults.fai_files, indexResults.dict_files, all_snps)
                filtered_indels = FilterIndels(reference_file, indexResults.fai_files, indexResults.dict_files, all_indels)

                if (params.annotation) {
                    annotated_snps = AnnotateSNPs(filtered_snps)
                    annotated_indels = AnnotateIndels(filtered_indels)
                }
    }

    // --- 5. PHYLOGENY ---
    if (params.phylogeny) {
        all_gvcfs = illumina_variants.gvcf_files.mix(ont_variants.gvcf_files).collect()
        snp_matrix_ch = Master_vcf(all_gvcfs, reference_file, indexResults.fai_files, indexResults.dict_files)

        if (params.annotation) {
            depth_files_collected = all_depth_files.collect()
            phylogeny_outputs = snp_matrix(snp_matrix_ch, depth_files_collected)
        } else {
            phylogeny_outputs = snp_matrix_no_annotate(snp_matrix_ch)
        }
    }
}

workflow.onComplete {
  println ( workflow.success ? "\nDone! Result files are in --> ./Outputs\n \
  If further analysis is required, bam alignments are in --> ./Outputs/bams\n \
  Phylogenetic tree and annotated merged variants are in --> ./Outputs/Phylogeny_and_annotation\n \
  Individual variant files are in --> ./Outputs/Variants/VCFs\n" \
  : "Oops .. something went wrong" )
}
