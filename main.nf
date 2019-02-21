#!/usr/bin/env nextflow


params.help = null


println """\

G E R M L I N E  V A R I A N T  D I S C O V E R Y  (S N V s  +  I N D E L S) 
============================================================================
Start    : $workflow.start

BWA      : BWA-0.7.15
Annovar  : Annovar-4.18
Samtools : Samtools-1.9
Strelka  : Strelka-2.9.7
Varscan  : Varscan-2.4.2
Picard   : Picard-2.18.15
GATK     : GenomeAnalysisTK-3.8-0

Results  : ${params.results}
"""
.stripIndent()


/* 
 * parse the input parameters.
 */

genome_fasta           = file(params.genome_fasta) 
genome_fasta_sa        = file(params.genome_fasta_sa)
genome_fasta_fai       = file(params.genome_fasta_fai)
genome_fasta_bwt       = file(params.genome_fasta_bwt)
genome_fasta_ann       = file(params.genome_fasta_ann)
genome_fasta_amb       = file(params.genome_fasta_amb)
genome_fasta_pac       = file(params.genome_fasta_pac)
genome_fasta_dict      = file(params.genome_fasta_dict)
mills_indel_file       = file(params.mills_indel_file)
known_indel_file       = file(params.known_indel_file)
known_dbsnps_file      = file(params.known_dbsnps_file)
known_dbsnps_1000_file = file(params.known_dbsnps_1000_file)


/*
 * read in the FASTQ files from s3 in pairs.
 */

raw_reads_ch = Channel.fromFilePairs("s3://noelnamai/data/${params.population}/*_{1,2}.fastq.gz")


/*
 * perform a variety of useful trimming tasks for Illumina data.
 */

process trimmomatic {

	tag "$sample"

	container "noelnamai/trimmomatic:0.38"

	cpus   = 2

	memory = "7.5 GB"

	when:
	sample == "ERR034520"

	input:
	set sample, file(reads) from raw_reads_ch

	output:
	set sample, file("${sample}.trimmed.paired.1.fq.gz"), file("${sample}.trimmed.paired.2.fq.gz") into trimmed_raw_reads_ch_1, trimmed_raw_reads_ch_2

	script:
	"""
	java -Xmx4g -jar /opt/trimmomatic-0.38.jar PE \
		-phred33 \
		${reads[0]} ${reads[1]} \
		${sample}.trimmed.paired.1.fq.gz ${sample}.trimmed.unpaired.1.fq.gz \
		${sample}.trimmed.paired.2.fq.gz ${sample}.trimmed.unpaired.2.fq.gz \
		ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
		LEADING:3 \
		TRAILING:3 \
		SLIDINGWINDOW:4:15 \
		MINLEN:36
	"""
} 


/*
 * perform a simple quality control check on raw sequence data coming from high throughput sequencing pipelines.
 */

process fastqc {

	tag "$sample"

	container "noelnamai/fastqc:0.11.3"

	cpus   = 2

	memory = "7.5 GB"

	publishDir "${params.results}/${params.population}/$sample", mode: "copy", overwrite: true

	input:
	set sample, file(forward), file(reverse) from trimmed_raw_reads_ch_1

	output:
	file("*.html") into fastqc_ch

	script:
	"""
	gunzip ${forward} ${reverse}

	cat "${forward.baseName}" "${reverse.baseName}" > ${sample}.fastq

	/usr/bin/fastqc ${sample}.fastq \
		--noextract \
		--threads ${task.cpus}

	mv ${sample}_fastqc.html ${sample}.fastqc.html
	"""
} 


/* 
 * read alignment via Burrows-Wheeler Aligner - MEM algorithm.
 */

process bwa_mem {

	tag "$sample"

	container "noelnamai/bwa:0.7.12"

	cpus   = 8

	memory = "30 GB"

	input:
	file genome_fasta     
	file genome_fasta_sa  
	file genome_fasta_fai 
	file genome_fasta_bwt 
	file genome_fasta_ann 
	file genome_fasta_amb 
	file genome_fasta_pac
	
	set sample, file(forward), file(reverse) from trimmed_raw_reads_ch_2

	output: 
	set sample, file("${sample}.sam") into bwa_aligned_sam_ch

	script:
	"""
	bwa mem -t ${task.cpus} -M ${genome_fasta} ${forward} ${reverse} > "${sample}.sam"
	"""
}


/* 
 * use Samtools to convert the SAM output from BWA to BAM format.
 */

process samtools_view {

	tag "$sample"

	container "noelnamai/samtools:1.9"

	cpus   = 2

	memory = "7.5 GB"

	input:
	set sample, file(sam_file) from bwa_aligned_sam_ch

	output: 
	set sample, file("${sam_file.baseName}.aligned.bam") into bwa_aligned_bam_ch

	script:
	"""
	samtools view -bS ${sam_file} > "${sam_file.baseName}.aligned.bam"
	"""
}


/* 
 * use Picard tools to add group to raw bam file.
 */

process picard_add_or_replace_read_groups {
	
	tag "$sample"

	container "noelnamai/picard:2.18.25"

	cpus   = 2

	memory = "7.5 GB"
	
	input:
	set sample, file(bam_file) from bwa_aligned_bam_ch
	
	output:
	set sample, file("${bam_file.baseName}.grouped.sorted.bam") into picard_added_group_bam_ch1, picard_added_group_bam_ch2
	
	script:
	"""
	java -Xmx4g -jar /opt/picard.jar AddOrReplaceReadGroups \
		I=${bam_file} \
		O="${bam_file.baseName}.grouped.sorted.bam" \
		SO=coordinate \
		RGLB=${sample} \
		RGPL=illumina \
		RGPU=${sample} \
		RGSM=${sample}
	"""
}


/* 
 * generate summary mapping statistics using Samtools.
 */

process samtools_flagstat {
	
	tag "$sample"

	container "noelnamai/samtools:1.9"

	cpus   = 2

	memory = "7.5 GB"

	publishDir "${params.results}/${params.population}/$sample", mode: "copy", overwrite: true
	
	input:
	set sample, file(bam_file_added_group) from picard_added_group_bam_ch1
	
	output:
	set sample, file("${sample}.mapping.statistic.txt") into samtools_flagstat_ch
	
	script:
	"""
	samtools flagstat ${bam_file_added_group} > "${sample}.mapping.statistic.txt"
	"""
}


/* 
 * use Picard tools to mark duplicates and create BAI from the BAM file.
 */

process picard_mark_duplicates {
	
	tag "$sample"

	container "noelnamai/picard:2.18.25"

	cpus   = 2

	memory = "7.5 GB"
	
	input:
	set sample, file(bam_file_added_group) from picard_added_group_bam_ch2
	
	output:
	set sample, file("${bam_file_added_group.baseName}.deduplicated.bam"), file("${bam_file_added_group.baseName}.deduplicated.bai") into picard_mark_duplicates_ch1, picard_mark_duplicates_ch2
	
	script:
	"""	
	java -Xmx4g -jar /opt/picard.jar MarkDuplicates \
		I=${bam_file_added_group} \
		O="${bam_file_added_group.baseName}.deduplicated.bam" \
		METRICS_FILE="${bam_file_added_group.baseName}.output.metrics" \
		VALIDATION_STRINGENCY=SILENT \
		CREATE_INDEX=true
	"""
}


/* 
 * GATK realigner target creator.
 */

process gatk_realigner_target_creator {

	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 8

	memory = "30 GB"

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	set sample, file(bam_file_marked_duplicate), file(bai_file_marked_duplicate) from picard_mark_duplicates_ch1
	
	output:
	set sample, file("${bam_file_marked_duplicate.baseName}.intervals.list") into gatk_target_creator_ch
	
	script:
	"""	
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T RealignerTargetCreator \
		-nt ${task.cpus} \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-o "${bam_file_marked_duplicate.baseName}.intervals.list" 		
	"""
}


/* 
 * GATK perform local realignment of reads around indels.
 */

process gatk_indel_realignment {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 8

	memory = "30 GB"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file mills_indel_file
	file known_indel_file
	
	set sample, file(interval_list) from gatk_target_creator_ch
	set sample, file(bam_file_marked_duplicate), file(bai_file_marked_duplicate) from picard_mark_duplicates_ch2
	
	output:
	set sample, file("${bam_file_marked_duplicate.baseName}.realigned.bam"), file("${bam_file_marked_duplicate.baseName}.realigned.bai") into gatk_indel_realignment_ch
	
	script:
	"""	
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T IndelRealigner \
		-R ${genome_fasta} \
		-I ${bam_file_marked_duplicate} \
		-known ${mills_indel_file} \
		-known ${known_indel_file} \
		-targetIntervals ${interval_list} \
		-o "${bam_file_marked_duplicate.baseName}.realigned.bam"
	"""
}


/* 
 * use GATK to detect systematic errors in base quality scores.
 */

process gatk_base_recalibrator {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 8

	memory = "30 GB"
	
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	file known_indel_file
	file mills_indel_file
	file known_dbsnps_file
	file known_dbsnps_1000_file	
	
	set sample, file(indel_realigned_bam), file(indel_realigned_bai) from gatk_indel_realignment_ch
	
	output:
	set sample, file(indel_realigned_bam), file(indel_realigned_bai), file("${indel_realigned_bam.baseName}.data.table") into gatk_base_recalibrator_ch
	
	script:
	"""	
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T BaseRecalibrator \
		-nct ${task.cpus} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		--knownSites ${known_dbsnps_1000_file} \
		--knownSites ${known_dbsnps_file} \
		--knownSites ${mills_indel_file} \
		--knownSites ${known_indel_file} \
		-o "${indel_realigned_bam.baseName}.data.table"
	"""
}


/* 
 * use GATK to write out sequence read data (for filtering, merging, subsetting etc).
 */

process gatk_print_reads {
	
	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 8

	memory = "30 GB"
		
	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	set sample, file(indel_realigned_bam), file(indel_realigned_bai), file(recalibration_table) from gatk_base_recalibrator_ch
	
	output:
	set sample, file("*.processed.bam"), file("*.processed.bai") into processed_bam_ch_1, processed_bam_ch_2
	
	script:
	"""
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T PrintReads \
		-nct ${task.cpus} \
		-R ${genome_fasta} \
		-I ${indel_realigned_bam} \
		-BQSR ${recalibration_table} \
		-o "${indel_realigned_bam.baseName}.processed.bam"
	"""
}


/* 
 * call Germline SNPs and indels via local re-assembly of Haplotypes.
 */

process haplotype_caller {

	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 2

	memory = "7.5 GB"

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	file known_dbsnps_file

	set sample, file(processed_bam), file(processed_bai) from processed_bam_ch_1

	output:
	set sample, file("${sample}.vcf") into haplotype_caller_ch

	script:
	"""
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T HaplotypeCaller \
		-R ${genome_fasta} \
		-I ${processed_bam} \
		--dbsnp ${known_dbsnps_file} \
		--out "${sample}.vcf"
	"""	
}


/* 
 * sort VCF files according to the order of the contigs in the header/sequence dictionary and then by coordinate. 
 */

process picard_sort_vcf {

	tag "$sample"

	container "noelnamai/picard:2.18.25"

	cpus   = 2

	memory = "7.5 GB"

	publishDir "${params.results}/${params.population}/$sample", mode: "copy", overwrite: true

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict
	
	set sample, file(original_vcf) from haplotype_caller_ch
	
	output:
	set sample, file("${original_vcf.baseName}.sorted.vcf") into picard_sort_vcf_ch

	script:
	"""
	java -Xmx4g -jar /opt/picard.jar SortVcf \
		I=${original_vcf} \
		O="${original_vcf.baseName}.sorted.vcf" \
		SEQUENCE_DICTIONARY=${genome_fasta_dict}
	"""
}


/* 
 * filter variant calls based on INFO and FORMAT annotations. 
 */

process gatk_variant_filtration {

	tag "$sample"

	container "broadinstitute/gatk3:3.8-0"

	cpus   = 2

	memory = "7.5 GB"

	publishDir "${params.results}/${params.population}/$sample", mode: "copy", overwrite: true

	input:
	file genome_fasta
	file genome_fasta_fai
	file genome_fasta_dict

	set sample, file(sorted_vcf) from picard_sort_vcf_ch

	output:
	set sample, file("${sorted_vcf.baseName}.filtered.vcf") into gatk_variant_filtration_ch

	script:
	"""
	java -Xmx4g -jar /usr/GenomeAnalysisTK.jar \
		-T VariantFiltration \
		-R ${genome_fasta} \
		-V ${sorted_vcf} \
		-window 35 \
		-cluster 3 \
		-filterName FS \
		-filter "FS > 30.0" \
		-filterName QD \
		-filter "QD < 2.0" \
		-o "${sorted_vcf.baseName}.filtered.vcf"
	"""
}


/* 
 * functionally annotate genetic variants detected from diverse genomes using ANNOVAR. 
 */

process annovar {

	tag "$sample"

	container "noelnamai/annovar:4.18"

	cpus   = 4

	memory = "15 GB"

	publishDir "${params.results}/${params.population}/$sample", mode: "copy", overwrite: true

	input:
	set sample, file(haplotype_caller_vcf) from gatk_variant_filtration_ch

	output:
	set sample, file("${sample}.annovar.hg19.multianno.txt") into annovar_ch

	script:
	"""
	perl /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene /opt/annovar/humandb/
	perl /opt/annovar/annotate_variation.pl -buildver hg19 -downdb cytoBand /opt/annovar/humandb/
	perl /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 /opt/annovar/humandb/
	perl /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 /opt/annovar/humandb/
	perl /opt/annovar/annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp30a /opt/annovar/humandb/

	perl /opt/annovar/table_annovar.pl ${haplotype_caller_vcf} /opt/annovar/humandb \
		-buildver hg19 \
		-out "${sample}.annovar" \
		-protocol refGene \
		-operation g \
		-nastring . \
		-vcfinput \
		--thread ${task.cpus} \
		--maxgenethread ${task.cpus} \
		-polish \
		-otherinfo

	mv ${sample}.annovar.hg19_multianno.txt ${sample}.annovar.hg19.multianno.txt
	"""
}
