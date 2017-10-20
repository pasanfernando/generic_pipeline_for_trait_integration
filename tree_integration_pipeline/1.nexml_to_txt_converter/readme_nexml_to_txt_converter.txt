This code converts the nexml character matrix downloaded from ontotrace into a tab delimited format, which is suitable to proceed through the pipeline. It correctly identifies the conflict states and writes ‘0 and 1’ for them. It also detects empty cell with missing data and writes ‘?’ for those cells. This code works for any matrix with any number of characters downloaded from ontotrace in nexml format.

input
the ontotrace matrix in nexml format (should be a .xml file)

outputs

tabdelemited_charactermatrix.txt  : The nexml matrix converted into tab delemited format. This is the input for the next code


