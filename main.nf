#!/usr/bin/env nextflow

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


Channel.fromFilePairs("s3://noelnamai/data/reads/*-WXS.read_{1,2}.fastq.gz")
	.map{file -> tuple(sampleName(file[0]), sampleType(file[0]), file[1])}
	.set{raw_reads_ch}


/* 
 * read alignment via Burrows-Wheeler Aligner - MEM algorithm and use of Samtools to convert the SAM output from BWA to BAM format.  
 */

process bwa_mem { 
	
	tag "$sampleName"

	container "biocontainers/bwa:0.7.15"

	input:
	file genome_fasta     
	file genome_fasta_sa  
	file genome_fasta_fai 
	file genome_fasta_bwt 
	file genome_fasta_ann 
	file genome_fasta_amb 
	file genome_fasta_pac
	
	set sampleName, sampleType, file(reads) from raw_reads_ch

	output: 
	set sampleName, sampleType, file("${sampleName}.${sampleType}.sam") into bwa_aligned_sam_wxs_ch

	script:
	"""
	bwa mem -t ${params.threads.mapping} -M ${genome_fasta} ${reads[0]} ${reads[1]} > "${sampleName}.${sampleType}.sam"
	"""
}