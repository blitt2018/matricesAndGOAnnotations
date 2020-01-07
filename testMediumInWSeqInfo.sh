#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=12:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --job-name="testMC05"
#SBATCH --mail-user=blitt@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module unuse /etc/modulefiles
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load blast-plus/2.6.0-py2-n35bbsh

myDB=/work/LAS/jernigan-lab/Ben/query_proteins/mouseAndHumanDB

myFastaFile=/work/LAS/jernigan-lab/Ben/query_proteins/mouseAndHumanMedium.faa

#create an empty directory and copy its path below
cd /work/LAS/jernigan-lab/Ben/testMediumInWSeqInfo

/work/LAS/jernigan-lab/Ben/blastMC05/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC05.out  -matrix BLOSUM90 
/work/LAS/jernigan-lab/Ben/blastMC10/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC10.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastMC13/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC13.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastMC20/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC20.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastMC22/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC22.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastTIP45/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out TIP45.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastBT1/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out BT1.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastD15/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out D15.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastMC26/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC26.out  -matrix BLOSUM90
/work/LAS/jernigan-lab/Ben/blastMC30/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out MC30.out  -matrix BLOSUM90

blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out BLOSUM62.out
