# Thesis' Results

In this directory are stored the diferent **results of Pablo Atienza' Bachelor Thesis**. The files could be separated into two groups, the CNV calling result tables, which end with '.tsv' and the SNP annotation results, which are the HTML files.

## CNV calling results

Each one of the six tables contains 39 thousand columns, one for each gene used for the CNV calling. The files have 14 rows, which are explained below. The files naming system is divided in three parts:

  1. **Final_table**: Indicates that the file has 14 rows in total, which indicate:
      * *Row 1*: Header with the genes names
      * *Rows 2-11*: One unnamed row for each individual in alphabetical order (HL0421, HL0438, HL0626, HL0672, HL0677, LL0409, LL0703, LL0720, LL0733, LL1010). In each position there can be three numbers:
          * 0: No CNV has been called
          * 1: A duplication has been called
          * -1: A deletion has been called
      * *Row 12*: Total number of individuals from both populations that have called a gene as a CNV.
      * *Row 13*: Total number of individuals from the Highlands population that have called a gene as a CNV.
      * *Row 14*: Total number of individuals from the Lowlands population that have called a gene as a CNV.
      
  2. **std/1kb**: The files under the tag 'std' have used the gene set for which the gene lenght has been left as original. The files under the tag '1kb' have the CNVs called using the lenght-modified gene set.
  3. **70Map/85Map/noMap**: Indication of the threshold for overlap with low mappability regions used for filtering the CNVs set. The '70Map' files are the most restrictive filtering, while 'noMap' indicates that no mappability filtering was applied.

## SNP annotation

The SNP annotation files are created using [SnpEff](http://pcingola.github.io/SnpEff/). They contain the results for the [Single Nucleotide Polymorphisms](https://github.com/TheRocinante-lab/PabloAtienza/blob/main/Results/Snp_summary.html) and for the [Indels](https://github.com/TheRocinante-lab/PabloAtienza/blob/main/Results/Indel_summary.html)

