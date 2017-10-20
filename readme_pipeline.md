This is a basic read me file for the pipeline. You can use the tutorial for step by step details on how to use the pipeline and visualize the character evolution in a phylogenetic tree. 

This pipeline transforms a character matrix obtained from Phenoscape Knowledge base (KB) into a version that can be efficiently integrated with an Open Tree phylogeny. 

	Optional: you can request for a meta-data file for the character your are interested in; this enables you to distinguish between inferred and asserted states (step 4 of the pipeline). Without the meta-data file you can still continue the pipeline, but you cannot differentially visualize the inferred vs asserted states. Refer to the tutorial for more details.

Please move the input files to the folder named ‘inputs’

To run the pipeline you need to use Python version 2 (2.7 or newer). This will not work on python 3. 

You need to install following external libraries as well
	1. networkx
	2. dendropy

If all the requirements are met, execute the main.py python script. The pipeline will be implemented and the output matrices can be found on the outputs folder.

tree_integration_pipeline folder contains the scripts for each step of the pipeline. You can look into the script and refer to the tutorial for more details.

1.‘preprocessed_tree.tre’: this is the pre-processed version of your input phylogeny from Open Tree of Life. This version must be used for the ancestral state reconstruction, because it does not have ott ids at the end of each species name. If you use the input version the ancestral state reconstruction will not run efficiently.

originaldatamatrix_taxalist.txt: contains a list of taxa in the original input data matrix

intermediate_matrices folder
