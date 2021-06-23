#!/usr/bin/env nextflow


/*
#==============================================
params
#==============================================
*/

params.saveMode = 'copy'
params.BAMfiles = "BAMfiles/*/*.bam"
params.Bedfile = "/home/patienza/Bedfiles/Mappability.bed"
params.Genes = "/home/patienza/RefGen/"
params.resultsDir = "CNVs_Map"

BAM_ch = Channel
        		.fromPath(params.BAMfiles)
                .map { file -> tuple(file.baseName, file) }


/*
#==============================================
Header log info
#==============================================
*/

log.info "========================================="
log.info "CNV calling with RD method"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="

/*
#==============================================
BAM files filtering for 30 MQ
#==============================================
*/

process Filtering {
	label 'medidum'

	input:
	tuple val(ID), path(BAM) from BAM_ch
	output:
	tuple val(ID), path("${ID}_30_filtered.bam") into BAM_ch

	script:

	"""
ml load SAMtools/1.9-GCC-8.2.0-2.31.1

samtools view -bq 30 ${BAM} > "${ID}_30_filtered.bam"
	"""
}

/*
#==============================================
BAM index
#==============================================
*/

process Indexing {
	label 'medidum'
	publishDir "FilteredBAM/$ID", mode: params.saveMode

	input:
	tuple val(ID), path(BAM) from BAM_ch
	output:
	tuple val(ID), path(BAM), path("${ID}.bai") into Index_ch

	script:

	"""
ml load SAMtools/1.9-GCC-8.2.0-2.31.1

samtools index ${BAM} > "${ID}.bai"
	"""

/*
#==============================================
Calculating genes mean with Bedtools
#==============================================
*/

process Genes_mean{
	label 'bigmem'
	
	input:
	tuple val(ID), path(BAM), path(Index) from Index_ch
	path Bedfile from params.Genes

	output:
	tuple val(ID), path(outFile) into Coverage_ch

	script:
	outFile = ID + '_GenesRD.txt'
	"""
	bedtools coverage -a $Bedfile -b $BAM -mean > $outFile   
	"""
}

/*
#==============================================
Comparison of means and CNV calling
#==============================================
*/

process Calling{
	publishDir "CNVs/$ID", mode: params.saveMode

	input:
	tuple val(ID), path(RDfile) from Coverage_ch

	output:
	tuple val(ID), path(CNVs) into CNV_ch

	script:
	CNVs = ID + '.cnv'
	"""
	ln -s /home/patienza/MexMaize/Scripts/Comparison.R Comparison.R
	Rscript Comparison.R $RDfile
	grep -v "Normal" $RDfile > $CNVs
	"""
}

/*
#==============================================
Comparison of means and CNV calling
#==============================================
*/

process Mappability {
	publishDir "Map_CNVs/$ID", mode: params.saveMode

	input:
	tuple val(ID), path(CNVs) from CNV_ch
	path Bedfile from params.Bedfile

	output:
	path outFile into Map_ch

	script:
	"""
	python /home/patienza/Programs/TCAG-WGS-CNV-workflow/compare_with_RLCR_definition.py $Bedfile "./" $CNVs 
	"""
}

process Analysis {
	publishDir "Map_CNVs/$ID", mode: params.saveMode

	input:
	path Map from Map_ch

	output:
	path outFile into Final

	script:
	"""
	Rscript Analysis.R $Map
	"""
}