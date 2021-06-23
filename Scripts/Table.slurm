#!/bin/bash

#SBATCH --output=Table-%J.%N.log # STDOUT file, %a for array, %A for job id, %N for hostname
#SBATCH --nodes=1 # request for cores
#SBATCH --mem 20G
#SBATCH --partition=fast # total run time limit in HH:MM:SS
#SBATCH --time=2:00:00 # total run time limit in HH:MM:SS
#SBATCH --mail-type=end # send email when job ends
#SBATCH --mail-user= p.atienza@upm.es


#srun Rscript Table_creator.R ~/RefGen/Genes/Names_B73v5.txt 

Var='_1kb_70Map'
n=0
for i in $(ls Map_CNVs/*/*_1kb_70_filtered.txt); do
  srun Rscript Table_modifier.R Genes_ID-base.txt $i $Var $n
  let n=$n+1
done

srun Rscript Final_row.R $Var
rm "Genes""$Var"".txt"

