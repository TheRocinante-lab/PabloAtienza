
# WORKFLOW FOR CNV CALLING ON MAIZE GENOME

The main scripts used for performning a CNV calling on maize genome based on a Read Depth method can be found. Inside the folder [ComplementaryScripts](https://github.com/TheRocinante-lab/PabloAtienza/tree/main/Scripts/ComplementaryScripts), the secondary scripts that are required for the workflow to work correctly are stored.

The script [Preprocessing.nf](https://github.com/TheRocinante-lab/PabloAtienza/blob/main/Scripts/Preprocessing.nf) is the first script to be run, since it contains all the steps from the trimming of the fasta files, to the creation of the final BAMfile.

Secondly, the [CNV calling](https://github.com/pabloati/Workflow_for_CNV_Calling/CNV_Calling) is performed by two methods:
  1. Extracting the mean read depth of each gene and comparing it against the total read depth of all the genes
  2. Using the tool [CNVnator](https://github.com/abyzovlab/CNVnator)
  
In order to produce a better and easiest way to interpret the results, the script [Table.sh](https://github.com/abyzovlab/CNVnator/Table.sh), creates a table with as many columns as genes are studied and as many rows as individuals. The genes that contain a duplication will be masked with a 1, deleted genes will contain a -1 and standard genes will contain a 0. The last row will be a sum of the number of indivuals that present structural variation (either duplication or deletion) for each gene studied.
