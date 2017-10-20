# this file works for the new format
#detects taxa with '0&1' and separates them according to taxa level
# prints the literature sources for each taxa : the output is intented to be reading only
# this code also replaces the 0&1 of internal nodes by '?'

input:
the xml files for fins: pectoral.xml or pelvic.xml (depending on which fin you want) 
missingremoved_matrix.txt: the original onto trace data matrix after the removal of missing taxa(without any processing)

Outputs:
intermediatecounts.txt: Contains the number of taxa that that has 0&1 states and their source paper where the data is coming from
conflicts_removed_datamatrix.txt: the data matrix where the 0&1 state of internal nodes are replaced by ‘?’. All the species level taxa still keeps 0&1 states
 
