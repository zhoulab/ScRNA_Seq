#!/bin/bash
#SBATCH --job-name=GSE130158_wDah_Female_ebi_download
#SBATCH --mail-type=ALL         
#SBATCH --mail-user=weihaotang@ufl.edu   
#SBATCH --mem-per-cpu=4gb
#SBATCH --cpus-per-task=2       
#SBATCH --qos=zhou-b               
#SBATCH --time=72:00:00            
#SBATCH --output=GSE130158_wDah_Female_ebi_download_%j.log   



wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/001/SRR8943131/SRR8943131.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/005/SRR8943125/SRR8943125.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/007/SRR8943127/SRR8943127.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/006/SRR8943126/SRR8943126.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/003/SRR8943133/SRR8943133.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/009/SRR8943119/SRR8943119.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/002/SRR8943132/SRR8943132.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/000/SRR8943120/SRR8943120.fastq.gz
wget -nc ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR894/001/SRR8943121/SRR8943121.fastq.gz
