# SPANDx
SPANDx - a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets

*****SPANDx now works with SLURM, SGE and PBS (Torque) resource managers (v2.7+) *****

SPANDx can also be run directly without a resource handler (set SCHEDULER=NONE in scheduler.config), although this is not recommended.


WHY USE SPANDX?

SPANDx is your one-stop tool for identifying SNP and indel variants in haploid genomes using NGS data. SPANDx performs alignment of raw NGS reads against your chosen reference genome or pan-genome, followed by accurate variant calling and annotation, and locus presence/absence determination. SPANDx produces handy SNP and indel matrices for downstream phylogenetic analyses. Annotated, genome-wide SNPs and indels are identified and output in human readable format. A presence/absence matrix is also generated to allow you to identify the core/accessory genome content across all your genomes.


USAGE: SPANDx.sh 
<parameters, required> 
-r <reference, without .fasta extension> 
[parameters, optional] 
-o [organism] 
-m [generate SNP matrix yes/no] 
-i [generate indel matrix yes/no] 
-a [include annotation yes/no] 
-v [Variant genome file. Name must match the SnpEff database] 
-s [Specify read prefix to run single strain] 
-t [Sequencing technology used Illumina/Illumina_old/454/PGM] 
-p [Pairing of reads PE/SE] 
-w [Window size in base pairs for BEDcoverage module]
-z [include tri- and tetra-allelic SNPs in the SNP matrix yes/no]

SPANDx by default expects reads to be paired end, Illumina data in the format: STRAIN_1_sequence.fastq.gz for the first pair and STRAIN_2_sequence.fastq.gz for the second pair. 
Reads not in this format will be ignored.
If your data are not paired, you must set the -p parameter to SE to denote unpaired reads. By default -p is set to PE.

SPANDx requires a reference file in FASTA format. 
For compatibility with all steps in SPANDx, FASTA files should conform to the specifications listed here: http://www.ncbi.nlm.nih.gov/BLAST/blastcgihelp.shtml
Note that the use of nucleotides other than A, C, G, or T is not supported by certain programs in SPANDx and these should not be used in reference FASTA files. 
In addition, Picard, GATK and SAMtools handle spaces within contig names differently. Therefore, please avoid using spaces or special characters (e.g. $/*) in contig names.

By default all reads in SPANDx format (i.e. strain_1_sequence.fastq.gz) in the current working directory will be processed. 
Sequence files will be aligned against the reference using BWA, alignments will be filtered and converted using SAMtools and Picard Tools.
SNPs and indels will be identified with GATK and coverage assessed with BEDtools.  SPANDx will then merge SNP and indel variants into matrices for phylogenetic resconstruction.
Identification of variants across multiple genomes involves several quality assessment steps in an attempt to minimize false positive variant calls 

Written by Derek Sarovich and Erin Price - Menzies School of Health Research, Darwin, Australia
Please send bug reports to mshr.bioinformatics@gmail.com or derek.sarovich@menzies.edu.au
If you find SPANDx useful and use it in published work, please cite! 
"SPANDx: a genomics pipeline for comparative analysis of large haploid whole genome re-sequencing datasets - BMC Research Notes 2014, 7:618"
