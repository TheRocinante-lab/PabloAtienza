#!/usr/bin/env nextflow

params.RefGen = "/home/patienza/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa"


process Split {
	label 'std'
	input:
	path Ref from params.RefGen

	output:
	path 'Splitted.fa' into Split_ch

	script:
	"""
	~/Programs/Seqbility/splitfa $Ref 35 > Splitted.fa 
	"""

}
process Splitting{
	
	label 'std'
	
	input:
	path Split from Split_ch

	output:
	path 'Chunk_*' into Splitted_ch

	script:
	"""
	split -l 218205270 Split > Chunk_
	"""
}
process Align {
	label 'long'

	input:
	path Ref from params.RefGen
	path Split from Splitted_ch

	output:
	tuple path(outFile), path(Split) into Sai_ch
	
	script:
	outFile = 'B73_v5_' + Split + '.sai'
	"""
 	module load BWA/0.7.17-GCC-8.2.0-2.31.1
	ln -s ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa.* ./
	bwa aln -R 1000000 -O 3 -E 3 $Ref $Split > $outFile
	"""
}


process Convert {
	label 'long'

	input:
	path Ref from params.RefGen
	tuple path(Sai), path(Split) from Sai_ch
	

	output:
	path outFile into Sam_ch

	script:
	outFile = 'B73_v5_' + Split + '.sam'
	"""
 	module load BWA/0.7.17-GCC-8.2.0-2.31.1
	ln -s ~/RefGen/B73_v5/Zm-B73-REFERENCE-NAM-5.0.fa.* ./
 	bwa samse $Ref $Sai $Split > $outFile
 	"""
}

Join = Sam_ch.collect()

process Merging {
	label 'long'

	input:
	path Files from Join

	output:
	path outFile into Merged

	script:
	outFile = 'B73_v5.sam'
	"""
	ml load SAMtools/1.9-GCC-8.2.0-2.31.1
	samtools merge $outFile $File
	"""
}

process Mask {
	label 'std'
	publishDir "~/RefGen/Map_DNA, mode: 'copy' "
	
	input:
	path SAM from Merged
	path Ref from params.RefGen

	output:
	path 'Zea_mays.AGPv4.dna_mm.toplevel.35kmer.fa' into Mask_ch

	script:
	"""
	~/Programs/Seqbility/gen_raw_mask.pl $SAM > rawMask_35.fa
	
	~/Programs/Seqbility/gen_mask -l 35 -r 0.5 rawMask_35.fa > mask_35_50.fa

	~/Programs/Seqbility/apply_mask_s mask_35_50.fa $Ref > Zea_mays.AGPv4.dna_mm.toplevel.35kmer.fa
	"""
}

