#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=01:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --job-name="blastP09936"
#SBATCH --mail-user=blitt@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module unuse /etc/modulefiles
module use /opt/rit/spack-modules/lmod/linux-rhel7-x86_64/Core/
module load blast-plus/2.6.0-py2-n35bbsh
module load parallel/20170322-36gxsog

myDB=/work/LAS/jernigan-lab/Ben/uniprot/human_proteome_db_container/human_proteome_db

myFastaFile=/work/LAS/jernigan-lab/Ben/P09936Experiment/P09936.faa 

#create an empty directory and copy its path below
cd /work/LAS/jernigan-lab/Ben/P09936Experiment
mkdir blastOutput
cd blastOutput

parallel -j16 "/work/LAS/jernigan-lab/Ben/blast{}/ncbi-blast-2.6.0+-src/c++/ReleaseMT/bin/blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt '6 qacc sacc evalue qstart qend sstart send qseq sseq' -out {}.out  -matrix BLOSUM90 " \
::: MC05 MC10 MC13 MC20 MC26 MC30 TIP45 D15 BT1

blastp  -num_threads 16  -db $myDB -query $myFastaFile -outfmt "6 qacc sacc evalue qstart qend sstart send qseq sseq" -out BLOSUM62.out

