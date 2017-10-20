This scripts detects the taxa with inferred presence and changes the state of the presence from ‘1’ to ‘2’. This new state is printed on new column. The original two columns for pelvic and pectoral fin will be kept without change.

input:	pectoralinferred.txt (the taxa list of inferred presence for pectoral fin)
	pelvicinferred.txt (the taxa list of inferred presence for pelvic fin)
	conflict_removed_datamatrix.txt (the input data matrix with 3 columns)
	
output: modified_inferredadded_matrix.txt ( the modified data matrix with two new 		columns for inferred presence state replacement)
	inferredstats: Prints the statistics for the number of inferred presence taxa that is transferred
	It will print out the taxa that was failed to be matched (due to name errors) find them and manually change the output matrix
