#!/usr/bin/env nextflow

params.help = null

println """\
S O M A T I C  S H O R T  V A R I A N T  D I S C O V E R Y  (S N V s  +  I N D E L S) 
=====================================================================================
Start      : $workflow.start

BWA        : BWA-0.7.15
Samtools   : Samtools-1.3.1
Picard     : Picard-2.3.0
GATK       : GenomeAnalysisTK-3.8-0
Strelka    : Strelka-2.9.7
Varscan    : VarScan-v2.4.2
"""
.stripIndent()


/* 
 * parse the input parameters.
 */

genome_fasta           = file(params.genome_fasta) 
genome_fasta_sa        = file(params.genome_fasta + ".sa")
genome_fasta_fai       = file(params.genome_fasta + ".fai")
genome_fasta_bwt       = file(params.genome_fasta + ".bwt")
genome_fasta_ann       = file(params.genome_fasta + ".ann")
genome_fasta_amb       = file(params.genome_fasta + ".amb")
genome_fasta_pac       = file(params.genome_fasta + ".pac")
genome_fasta_dict      = file("${genome_fasta.parent}/ucsc.hg19.dict")
cosmic_file            = file(params.cosmic_file)
mills_indel_file       = file(params.mills_indel_file)
known_indel_file       = file(params.known_indel_file)
known_dbsnps_file      = file(params.known_dbsnps_file)
known_dbsnps_1000_file = file(params.known_dbsnps_1000_file)


/* helper functions, given a file path returns the file name region */

def sampleName(file) {

	def tokens  = file.tokenize("-")
	String name = tokens[0]

	return name
}

def sampleType(file) {

	def tokens  = file.tokenize("-")
	String type = tokens[1]

	switch(type) {

		case "T":
			return "tumor"
			break;
		case "N":
			return "normal"
			break;
	}
}


// raw_reads_ch = Channel.fromFilePairs("s3://noelnamai/data/reads/*-WXS.read_{1,2}.fastq.gz")
// 	.map{file -> tuple(sampleName(file[0]), sampleType(file[0]), file[1])}


raw_reads_ch = Channel.from( 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)


/* 
 * read alignment via Burrows-Wheeler Aligner - MEM algorithm and use of Samtools to convert the SAM output from BWA to BAM format.  
 */

process bwa_mem { 
	
	tag "$number"

	container "biocontainers/bwa:0.7.15"

	publishDir "${params.results}", mode: "move", overwrite: true

	input:
	file genome_fasta
	file genome_fasta_sa
	file genome_fasta_fai
	file genome_fasta_bwt
	file genome_fasta_ann
	file genome_fasta_amb
	file genome_fasta_pac
	
	val number from raw_reads_ch

	output: 
	file("${number}.txt") into bwa_aligned_sam_wxs_ch

	script:
	"""
	echo "$number" > "${number}.txt"
	"""
}