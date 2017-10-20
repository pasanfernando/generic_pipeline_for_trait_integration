The OTnaming_pipeline.py script is used to convert the VTO data matrix into open tree version.
It will use combined approach that includes name conversion using NCBI ids and direct name matching.
The script can also generate statistics for the final conversion.

Inputs

finalVTOmatrix.txt: this is the final propagated matrix and has VTO taxa
all_tips.txt: A list of all tips (species from final open tree file)
taxonomy.tsv: The data file downloaded from open tree database (contains all the information about all the taxa in open tree. Used to extract NCBI ids)
vtonewfinal.owl: The VTO ontology file

Outputs
finalfullopentree_matrix.txt: The final matrix. Contains open tree taxa. Ready for the mapping. This file contains all the open tree taxa that is in the tree file. Even the ones without data.

finalopentree_matrix_onlydata.txt: This final matrix contains only the taxa with data. It can be used for the mapping as well

finalmismatchedlist_andstats.txt: This file contains all the statistics about matching and mismatching taxa. Also prints out the mismatching taxa and counts how many extinct and sp. (inaccurate) names are there in the final mismatched taxa list.
Also prints the mismatched taxa with finloss separately

finalmatchedvtlist.txt:list of all matched taxa in VTO form
