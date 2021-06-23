#!/usr/bin/env nextflow


/*
#==============================================
params
#==============================================
*/

params.saveMode = 'copy'
params.filespath = "FastaFiles/*/*{1,2}.fastq"
params.resultsDir = "BAMfiles"
params.refGen = "/home/patienza/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa"

ch_trim = Channel.fromFilePairs(params.filespath)

/*
#==============================================
Header log info
#==============================================
*/

log.info "========================================="
log.info "GATK Best Practices for CNV Preprocessing"
log.info "Nextflow Version:	$workflow.nextflow.version"
log.info "Command Line:		$workflow.commandLine"
log.info "========================================="


/*
#==============================================
trimmomatic	
#==============================================
*/


process Trimming {
	tag "$ID"
	label "long"
	publishDir "Trimmed/$ID", mode: params.saveMode
	
	input:
	tuple val(genName), file(reads) from ch_trim

 	output:
	tuple val(genName), path(fq_1_paired), path(fq_2_paired) into Paired_trim
	tuple val(genName), path(fq_1_unpaired), path(fq_2_unpaired) into Unpaired_trim
	
	script:
	ID = genName.take(6)
	fq_1_paired = genName + '_R1_paired.fastq'
   	fq_1_unpaired = genName + '_R1_unpaired.fastq'
    	fq_2_paired = genName + '_R2_paired.fastq'
    	fq_2_unpaired = genName  + '_R2_unpaired.fastq'
	"""
	ml load Trimmomatic/0.38-Java-1.8
	java -jar \$EBROOTTRIMMOMATIC/trimmomatic-0.38.jar \
	PE \
	-threads 10 \
	-phred33 \
	$reads \
	$fq_1_paired \
	$fq_1_unpaired \
	$fq_2_paired \
	$fq_2_unpaired \
	LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 
	"""
}


/*
#==============================================
BWA mem alignment
#==============================================
*/


process PairedEnd {
	label 'long'
	tag "$ID"

	input:
	tuple val(genName), file(reads) from Paired_trim
	path Ref from params.refGen

 	output:
	tuple val(genName), path(outFile) into Aligned_paired
	
	script:
	ID = genName.take(6)
	outFile = genName + '_B73v5_PE.bam'
	Header = '@RG\\tID:NRGENE\\tSM:' + ID + '\\tPL:Illumina\\tLB:no\\tPU:unit1'
	"""
 	module load BWA/0.7.17-GCC-8.2.0-2.31.1
	ml load SAMtools/1.9-GCC-8.2.0-2.31.1
 	ln -s ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa.* ./
bwa mem -t 10 -M -R "$Header" $Ref $reads | samtools view -bS - > $outFile
	"""
}


process UnpairedReads {
	label 'long'
	tag "$ID"

	input:
	tuple val(genName), file(reads) from Unpaired_trim
	path Ref from params.refGen

 	output:
	tuple val(genName), path("$outFile1"), path("$outFile2") into Aligned_unpaired
	
	script:
	ID = genName.take(6)
	outFile1 = genName + '_R1_B73v5_SE.bam'
	outFile2 = genName + '_R2_B73v5_SE.bam'
	Header = '@RG\\tID:NRGENE\\tSM:' + ID + '\\tPL:Illumina\\tLB:no\\tPU:unit1'

	"""
 	module load BWA/0.7.17-GCC-8.2.0-2.31.1
 	module load SAMtools/1.9-GCC-8.2.0-2.31.1
 	ln -s ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa.* ./
bwa mem -t 10 -M -R "$Header" $Ref ${reads[0]} | samtools view -bS - > $outFile1
bwa mem -t 10 -M -R "$Header" $Ref ${reads[1]} | samtools view -bS - > $outFile2
	"""
}
/*
#==============================================
# Picard-TOOLS sorting
#==============================================
*/

process Sort_paired {
	label 'medium'
	tag "$ID"

	input:
	tuple val(genName), path(pairedFile) from Aligned_paired
	

	output:
	tuple val(ID), path(outFile) into Sorted_paired

	script:
	ID = genName.take(6)
	outFile = genName + '_PE_sorted.bam'
	"""
	java -Xmx20g -XX:ParallelGCThreads=10 -jar ~/Programs/PicardTools/picard.jar SortSam \
	MAX_RECORDS_IN_RAM=2000000 \
	INPUT=$pairedFile \
	OUTPUT=$outFile \
	SORT_ORDER=coordinate
	"""
}

process Sort_unpaired {
	label 'medium'
	tag "$ID"

	input:
	tuple val(genName), path(File1), path(File2) from Aligned_unpaired

	output: 
	tuple val(ID), path(outFile1), path(outFile2) into Sorted_unpaired

	script:
	ID = genName.take(6)
	outFile1 = genName + '_SE_sorted_1.bam'
	outFile2 = genName + '_SE_sorted_2.bam'
	"""
	java -Xmx20g -XX:ParallelGCThreads=10 -jar ~/Programs/PicardTools/picard.jar SortSam \
	MAX_RECORDS_IN_RAM=2000000 \
	INPUT=$File1 OUTPUT=$outFile1 \
	SORT_ORDER=coordinate

	java -Xmx20g -XX:ParallelGCThreads=10 -jar ~/Programs/PicardTools/picard.jar SortSam \
	MAX_RECORDS_IN_RAM=2000000 \
	INPUT=$File2 OUTPUT=$outFile2 \
	SORT_ORDER=coordinate
	"""
}

/*
#==============================================
# Merging data using SAMtools
#==============================================
*/

Paired_grouped = Sorted_paired.groupTuple()
Unpaired_grouped = Sorted_unpaired.groupTuple()
Sorted_joined = Paired_grouped.join(Unpaired_grouped)

process Merging {
	label 'medium'
	tag "$ID"

	input:
	tuple val(genName), path(pairedFile), path(unpairedFile1), path(unpairedFile2) from Sorted_joined

	output:
	tuple val(ID), path(outFile) into Merged

	script:
	ID = genName.take(6)
	outFile = ID + '_merged.bam'
	"""
	ml load SAMtools/1.9-GCC-8.2.0-2.31.1
	samtools merge $outFile $pairedFile $unpairedFile1 $unpairedFile2 
	"""
}

/*
#==============================================
# Removing duplicates using Picard Tools
#==============================================
*/

process Rem_Dup {
	publishDir "Raw_BAMfiles/$ID", mode: params.saveMode, pattern: "*_metrics.txt"
	publishDir "Raw_BAMfiles/$ID", mode: params.saveMode, pattern: "*_dedup.bam"
	label 'bigmem'
	tag "$ID"

	input: 
	tuple val(ID), path(BAMFile) from Merged

	output:
	tuple val(ID), path(outFile), path(metricFile) into Deduped

	script:
	outFile = ID +'_dedup.bam'
	metricFile = ID + '_metrics.txt'
	"""
 	java -Xmx32g -XX:ParallelGCThreads=10 -Djava.io.tmpdir=/scratch/ -jar ~/Programs/PicardTools/picard.jar \
	MarkDuplicates \
	INPUT=$BAMFile OUTPUT=$outFile \
	MAX_RECORDS_IN_RAM=2000000 \
	REMOVE_DUPLICATES=true \
	METRICS_FILE=$metricFile 
 	"""
}

/*
#==============================================
# Indexing using Picard Tools
#==============================================
*/

process Indexing {
	label 'bigmem'
	tag "$ID"

	input:
	tuple val(ID), path(inFile), path(metricFile) from Deduped

	output:
	tuple val(ID), path(inFile), path(indexFile) into Indexed

	script:
	indexFile = ID + '_dedup.bai'
	"""
	java -Xmx200g -XX:ParallelGCThreads=10 -jar -Djava.io.tmpdir=/scratch/ ~/Programs/PicardTools/picard.jar \
	BuildBamIndex \
	INPUT=$inFile OUTPUT=$indexFile
	"""
}

/*
#==============================================
# Indel Realignment using GATK
#==============================================
*/
process Target_creator{
	label 'bigmem'
	tag "$ID"

	input:
	tuple val(ID), path(BAMfile), path(indexFile) from Indexed
	path Ref from params.refGen

	output:
	tuple val(ID), path(BAMfile), path(indexFile), path(outFile) into Targeted

	script:
	outFile = ID + '-forIndelReaigner.intervals'
	"""
	module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
	ln -sfn ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.* ./
	java -Xmx200g -jar ~/Programs/GATK3.8/GenomeAnalysisTK.jar -T \
	RealignerTargetCreator \
	-R $Ref \
	-I $BAMfile -o $outFile
	"""
}

process Indel_Realigner {
	label 'bigmem'
	tag "$ID"
	publishDir "$params.resultsDir/$ID/", mode: params.saveMode

	input:
	tuple val(ID), path(BAMfile), path(indexFile), path(targetsFile) from Targeted
	path Ref from params.refGen

	output:
	path outFile into Indel_Realigned

	script:
	outFile = ID + '_indelrealigned.bam'
	"""
	module load GATK/4.1.2.0-GCCcore-8.2.0-Java-1.8
	ln -sfn ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.* ./
	java -Xmx200g -jar ~/Programs/GATK3.8/GenomeAnalysisTK.jar -T \
	IndelRealigner \
	-R $Ref \
	-I $BAMfile  -targetIntervals $targetsFile \
	-o $outFile
	"""
}

