This scripts propagates data of internal nodes to their species. If the species already have data for a specific internal node that data will be kept. All the other species without data will be added to data file. For now, the propagation only considers internal taxa up to family level.

input:	modified_inferredadded_matrix.txt (the input data matrix with 5 columns)
	vtonewfinal.owl(VTO ontology data file; the syntax for one species was 		problematic in the direct download. This was changed in this file.
	
output: finalVTOmatrix.txt (the final matrix with 7 columns:two new columns were added to indicate wether the species is propagated or not; ‘1’:propagated, 0: not propagated)
	propagationstatistics.txt (gives statistics about the original data file, propagation step and the new data matrix. It also gives the counts of propagated taxa for families and genera separately)
	finalVTOspecieslist.txt (list of species in the finalVTOmatrix. can be used for comparison purposes)
	famlilyand_genera_propahated.txt  statistics for families and genera that propagated data
