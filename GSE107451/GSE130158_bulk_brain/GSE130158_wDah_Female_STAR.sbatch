#!/bin/bash
#SBATCH --job-name=GSE130158_wDah_Female_STAR
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=weihaotang@ufl.edu   
#SBATCH --mem-per-cpu=8gb
#SBATCH --cpus-per-task=6        
#SBATCH --qos=zhou-b               
#SBATCH --time=72:00:00            
#SBATCH --output=GSE130158_wDah_Female_STAR_%j.log   

module load star 
#build star genome index
#STAR --runThreadN 23 --runMode genomeGenerate --genomeDir ./STARgenome643 --genomeFastaFiles ./dmel-all-chromosome-r6.43.fasta --sjdbGTFfile ./dmel-all-r6.43.gtf


cd /orange/zhou/projects/aging/GSE130158_Brain
for i in *.fastq.gz;do
names=${i%.fastq*} 
STAR --readFilesCommand zcat --genomeDir /orange/zhou/projects/aging/SDIII/CuffLinks/STARgenome643 --sjdbGTFfile /orange/zhou/projects/aging/SDIII/CuffLinks/dmel-all-r6.43.gtf \
        --runThreadN 24 --runMode alignReads --outSAMtype BAM SortedByCoordinate --outBAMsortingBinsN 200 --limitBAMsortRAM 16013050982\
        --readFilesIn ./${names}.fastq.gz \
         --outFileNamePrefix /orange/zhou/projects/aging/GSE130158_Brain/${names}
done





