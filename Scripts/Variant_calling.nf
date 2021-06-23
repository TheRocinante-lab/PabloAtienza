#!/usr/bin/env nextflow


/*
#==============================================
params
#==============================================
*/

params.saveMode = 'copy'
params.bamFiles = "Raw_BAMfiles/*/*_dedup.bam"
params.refGen = "/home/patienza/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa"
params.chrFile = ""
BAM_ch = Channel.fromPath(params.bamFiles)
                 .map { file -> tuple(file.baseName.take(6), file) }

BAM_ch.into { Test_ch; Filter_ch }
/*
#==============================================
Header log info
#==============================================
*/

log.info "========================================="
log.info "GATK Best Practices for Variant Calling"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="


process Ini_Var_Hap_Caller {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(bam) from Test_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(vcf) into Raw_variants_snp_ch, Raw_variants_indel_ch

	script:

	vcf = ID + "_raw_variants.vcf"

	"""
	gatk HaplotypeCaller \
	-R $ref \
	-I $bam \
	-O $vcf
	"""

}

process Ini_var_SelectSNP {
	tag "$ID"
	label 'fast'
	
	input:
	tuple val(ID), path(vcf) from Raw_variants_snp_ch
	path ref from params.refGen


	output:
	tuple val(ID), path(snp_out) into Raw_snp_ch

	script:
	snp_out = ID + "_raw_snps.vcf"
	"""

	gatk SelectVariants \
	-R $ref \
	-V $vcf \
	-select-type SNP \
	-O $snp_out
	"""
}

process Ini_var_SelectINDEL {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(vcf) from Raw_variants_indel_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(indel_out) into Raw_indel_ch

	script:

	indel_out = ID + "_raw_indels.vcf"
	"""

	gatk SelectVariants \
	-R $ref \
	-V $vcf \
	-select-type INDEL \
	-O $indel_out	
	"""
}

process Ini_var_FilterSNP {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(snp) from Raw_snp_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(out) into IV_Filtered_snp_ch

	script:
	out = ID + "_filtered_snps.vcf"
	"""

	gatk VariantFiltration \
	-R $ref \
	-V $snp \
	-O $out \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum"
	"""
}

process Ini_var_FilterINDEL {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(indel) from Raw_indel_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(out) into IV_Filtered_indel_ch

	script:
	out = ID + "_filtered_indels.vcf"
	"""

	gatk VariantFiltration \
	-R $ref \
	-V $indel \
	-O $out \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	"""
}

process Ini_var_FinalSNP {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(snp) from IV_Filtered_snp_ch

	output:
	tuple val(ID), path(out) into BQSR_snp_ch

	script:
	out = ID + "_bqsr_snps.vcf"
	"""

	gatk SelectVariants \
	--exclude-filtered \
	-V $snp \
	-O $out
	"""

}

process Ini_var_FinalINDEL {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(indel) from IV_Filtered_indel_ch

	output:
	tuple val(ID), path(out) into BQSR_indel_ch

	script:
	out = ID + "_bqsr_indels.vcf"
	"""

	gatk SelectVariants \
	--exclude-filtered \
	-V $indel \
	-O $out
	"""
}

BQSR_ch = BQSR_snp_ch.join(BQSR_indel_ch)
Unrecal_ch = Filter_ch.join(BQSR_ch)

process Before_Recal_Table {
	tag "$ID"

	label 'fast'

	input:
	tuple (ID), path(bam), path(bqsr_snp), path(bqsr_indel) from Unrecal_ch
	path ref from params.refGen


	output:
	tuple val(ID), path(bam), path(bqsr_snp), path(bqsr_indel), path(table) into For_recalling_ch
	
	script:
	table = ID + "_recal_data.table"
	"""

	gatk BaseRecalibrator \
	-R $ref \
	-I $bam \
	--known-sites $bqsr_snp \
	--known-sites $bqsr_indel \
	-O $table
	"""
}

process Apply_Base {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(bam), path(bqsr_snp), path(bqsr_indel), path(table) from For_recalling_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(recal_bam), path(bqsr_snp), path(bqsr_indel) into Recal_ch
	tuple val(ID), path(table) into Before_table_ch
	tuple val(ID), path(recal_bam) into BAM_recalled_ch
	
	script:
	recal_bam = ID * "_recal.bam"
	"""

	gatk ApplyBQSR \
	-R $ref \
	-I $bam \
	-bqsr $table \
	-O $recal_bam
	"""
}

process After_Recal_Table {
	tag "$ID"
	label 'fast'

	input:
	tuple val(ID), path(recal_bam), path(bqsr_snp), path(bqsr_indel) from Recal_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(table) into After_table_ch

	script:
	table = ID + "_post_recal_data.table"
	"""

	gatk BaseRecalibrator \
	-R $ref \
	-I $bam \
	--known-sites $bqsr_snp \
	--known-sites $bqsr_indel \
	-O $table
	"""
}

Tables_ch = Before_table_ch.join(After_table_ch)

process Analisis_Recaling {
	tag "$ID"
	publishDir "BSQR/", mode: params.saveMode, pattern: "*.pdf"
	label 'fast'

	input:
	tuple val(ID), path(table_before), path(table_after) from Tables_ch

	output:
	path PDF into PDF_ch

	script:
	PDF = ID + "_recalibration_plots.pdf"
	"""
	gatk AnalyzeCovariates \
	-before $table_before \
	-after $table_after \
	-plots $PDF
	"""
}

/*
#==============================================
# Call variant per sample
#==============================================
*/

process Variant_Caller {
	tag "$ID"
	publishDir "GVCF", mode: params.saveMode
	label 'fast'

	input:
	tuple val(ID), path(bam) from BAM_recalled_ch
	path ref from params.refGen

	output:
	tuple val(ID), path(out) into Vcf_ch

	script:
	out = ID + '.g.vcf.gz'
	"""

	gatk --java-options "-Xmx4g" HaplotypeCaller \
	-R $ref \
	-I $bam \
	-O $out\
	-ERC GVCF
	"""
}

List_ch = Vcf_ch.collectFile(name:"maize.sample_map.txt", newLine: true, storeDir: "~/MexMaize/GVCF/"){it[0] + '\n' + it[1]}

process Assembly {
	tag "All"
	publishDir "Maize_db", mode: params.saveMode
	label 'fast'

	input:
	path list from List_ch
	path interval from params.chrFile

	output:
	val db_dir into DB_ch

	script:
	db_dir = "/home/patienza/MexMaize/Maize_db"
	"""

	gatk --java-options "-Xmx8g -Xms4g" GenomicsDBImport \
	--sample-name-map $list \
	--genomicsdb-workspace-path $db_dir \
	-L $interval 
	"""
}

Chr_ch = Channel.from(1..10)

process Unfiltered_vars {
	tag "Chr $num"
	publishDir "Variants/Unfiltered", mode: params.saveMode
	label 'fast'

	input:
	val num from (1..10)
	each val db_dir from DB_ch
	path ref from params.refGen

	output:
	tuple val(num), path(out) into Unfiltered_geno_ch

	script:
	out = 'geno_' + $num + '_unfiltered.vcf.gz'
	"""

	gatk --java-options "-Xmx4g" GenotypeGVCFs \
	-R $ref \
	-V gendb://$dir \
	-L $num \
	-O $out
	"""
}

/*
#==============================================
# Filtering
#==============================================
*/

Unfiltered_geno_ch.into {Unfiltered_geno_snp_ch; Unfiltered_geno_indel_ch}

process Selecting_SNPs{
	tag"Chr $num"
	label 'fast'

	input:
	tuple val(num), path(vars) from Unfiltered_geno_snp_ch
	path $ref from params.refGen

	output:
	tuple val(num), path(out) into Unfiltered_snp_ch

	script:
	out = 'geno_' + num + '_unfiltered_snps.vcf'
	"""

	gatk SelectVariants \
	-R $ref \
	-V $vars \
	-select-type SNP \
	-O $out \
	"""
}

process Selecting_INDELSs{
	tag "Chr $num"
	label 'fast'

	input:
	tuple val(num), path(vars) from Unfiltered_geno_indel_ch
	path $ref from params.refGen

	output:
	tuple val(num), path(out) into Unfiltered_indel_ch

	script:
	out = 'geno_' + num + '_unfiltered_indel.vcf'
	"""
	
	gatk SelectVariants \
	-R $ref \
	-V $vars \
	-select-type INDEL \
	-O $out \
	"""
}

process Filter_SNPs {
	tag "Chr $num"
	label 'fast'

	input:
	tuple val(num), path(snps) from Unfiltered_snp_ch
	path ref from params.refGen

	output:
	tuple val(num), path(out) into Filtered_snp_ch

	script:
	out = 'geno_' + num + '_filtered_snps.vcf'
	"""

	gatk VariantFiltration \
	-R $ref \
	-V $snps \
	-O $out \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"
	"""
}

process Filter_INDELs {
	tag "Chr $num"
	label 'fast'

	input:
	tuple val(num), path(indels) from Unfiltered_indel_ch
	path ref from params.refGen

	output:
	tuple val(num), path(out) into Filtered_indel_ch

	script:
	out = 'geno_' + num + '_filtered_indels.vcf'
	"""

	gatk VariantFiltration \
	-R $ref \
	-V $indels \
	-O $out \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" 
	"""
}

process Final_SNPs{
	tag "Chr $num"
	publishDir "Variants/Filtered", mode: params.saveMode
	label 'fast'


	input:
	tuple val(num), path(snps) from Filtered_snp_ch

	output:
	path out into Final_chr_snp_ch

	script:
	out= 'geno_' + num + '_final_snp.vcf'
	"""
	gatk SelectVariants \
	--exclude-filtered \
	-V $snps \
	-O $out
	"""
}

process Final_INDELs{
	tag "Chr $num"
	publishDir "Variants/Filtered", mode: params.saveMode
	label 'fast'

	input:
	tuple val(num), path(indels) from Filtered_indel_ch

	output:
	path out into Final_chr_indel_ch

	script:
	out= 'geno_' + num + '_final_indel.vcf'
	"""
	gatk SelectVariants \
	--exclude-filtered \
	-V $indels \
	-O $out
	"""
}


