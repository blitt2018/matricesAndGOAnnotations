# BLAST Matrices and Gene Ontology (GO) Annotations 
This project is focused on using GO annotations as a means to confirm homology between proteins. In particular, new proteins were predicted using custom scoring matrices, and the goal was to see if these new proteins were homologous using their GO annotations. The pipeline is detailed below. 

## BLASTing 
First, it was necessary to BLAST proteins using the novel matrices. Code is not provided for this; however, BLAST had to be recompiled for every different matrix used. After this recompilation, bash scripts were used to BLAST identical queries against the same database. This was done on an HPC cluster using BLAST 2.6. Example found in testMediumInWSeqInfo.sh. 

## Visualizing BLAST Output Quantity 
One simple but crucial peice of information is the amount of hits returned by versions of BLAST employing different scoring matrices. This is the start of classifying what a "better" BLAST would look like. One consideration is that BLAST outputs multiple hits for the same protein, if the protein is matched at various points of its sequence. Visualizations were made that both accounted for, and did not account for this factor. Example output shown below: 

![Amount of hits w/o doubles](https://github.com/blitt2018/matricesAndGOAnnotations/blob/master/blastOutAmpGroupsNoDoubles.png)
![Amount of hits w/ doubles](https://github.com/blitt2018/matricesAndGOAnnotations/blob/master/blastOutAmpGroupsWDoubles.png)

## Retreiving GO Annotations 
After generating BLAST output from different scoring matrices, the next step was to retreive GO annotations for each output protein. This was done along two pathways, depending on the source of proteins used when BLASTing. 

- If using the NCBI RefSeq database, the easiest way to get GO Annotations was via progrommatic queries to the Uniprot API. For RefSeq IDs, GOsToFileHierarchy.py was used. 

- If using UniProt IDs, Gene Ontology information can be found using a crossreference file from the UniProt website. The script, GOsToFileHierarchyLocalParallel.py, uses the crossreference file and multithreading to annotate proteins with their GO functions. 

These files both create a hierarchy where each leaf is a file containing all of the GO annotations for one BLAST query used at the start of the pipeline.

## Plotting GO Term Frequencies 
One useful method of comparing BLOSUM62 with our novel scoring matrices is to plot the frequencies of each GO term produced by either of the matrices. This shows how often the novel scoring matrix in question missed a GO term produced by BLOSUM62, how often it found the same GO term as BLOSUM62, and how often it found a GO term that BLOSUM62 did not. One important note is that when processing BLAST output, multiple hits with the same protein ID were combined and considered one hit. 

To generate this figure, the protein YP_003024038 was input to BLAST and the output was collected for different scoring matrices. The GO annotations for the BLAST output were collected and plotted for each matrix. These graphs are interactive and can be zoomed in on to inspect various sections. 
![GO frequency plots](https://github.com/blitt2018/matricesAndGOAnnotations/blob/master/reproducePipelineMediumInputAttempt/GOFreqsYP_003024038ZoomedIn.png)

Here it is shown that some matrices have a unique set of GO terms that appear more often than the rest. Ideally, these peaks would be from proteins with functions that are more specific than the B62 output, but still closely related enough to be reasonable outputs. For this reason, the script used to generate the visual above also has command line output. This output shows whether any GO terms in these peaks are children of functions annotated to proteins given by BLOSUM62. 

*Usage*
- Use the script plotAllQueriesIndividualQueryPipeline.py to produce graphs for every query in the blast output directory. 
- Uncomment line 144 to automatically zoom in on the relevant part of each graph (won't all graphs on the same scale).
- Adjust TERM_LIMIT to control the amount of GO terms unique to the novel matrix shown. 
- Adjust LEVELS_TO_SEARCH to control how many levels are searched in the GO hierarchy for parents output by BLOSUM62. 
