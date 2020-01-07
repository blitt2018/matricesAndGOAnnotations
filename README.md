# BLAST Matrices and Gene Ontology (GO) Annotations 
This project is focused on using GO annotations as a means to confirm homology between proteins. In particular, new proteins were predicted using custom scoring matrices, and the goal was to see if these new proteins were homologous using their GO annotations. The pipeline is detailed below. 

## BLASTing 
First, it was necessary to BLAST proteins using the novel matrices. Code is not provided for this; however, BLAST had to be recompiled for every different matrix used. After this recompilation, bash scripts were used to BLAST identical queries against the same database. This was done on an HPC cluster using BLAST 2.6. Example found in testMediumInWSeqInfo.sh. 

## Retreiving GO Annotations 
After generating BLAST output from different scoring matrices, the next step was to retreive GO annotations for each output protein. This was done along two pathways, depending on the source of proteins used when BLASTing. 

- If using the NCBI RefSeq database, the easiest way to get GO Annotations was via progrommatic queries to the Uniprot API. For RefSeq IDs, GOsToFileHierarchy.py was used. 

- README for Uniprot IDs needs to be written 

