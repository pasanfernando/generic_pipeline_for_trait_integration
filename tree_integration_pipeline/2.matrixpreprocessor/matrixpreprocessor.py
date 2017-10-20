# author Pasan Fernando
# Date 06/13/2016
# Takes the tab delemited matrix that is converted from nexml format.
# prints all the taxa into a separate file
# Also removes the missing taxa ( ones with two ? for both pelvic and pectoral fins) and prints the statistics
__author__ = 'pasan fernando'
######################################################################################################################
import sys
import re
import networkx as nx
import collections
import os

# all the paths should be in relation to the main.py script
sys.path.insert(0, os.path.join("tree_integration_pipeline","lib"))
import methods as mt

###### Reading the Vertebrate Taxonomy ontology (VTO) file to store relationships

name,namer,rank,G =mt.VTOreader()

######################################################################################################################


missinglist =[]
alltaxa =[]
remainingspecies ={}

######Read the tab delimited data matrix
#in1 = raw_input('enter the input matrix:')

d = open(os.path.join("intermediate_matrices","tabdelemited_charactermatrix.txt"), 'r')

for line in d:
    if (line != '\n') and ('taxa_name' not in line):
        line = line.strip('\n')
        a = line.split('\t')
        a[0] = a[0].strip('\'')
        # following replacements removes the conventional naming errors in VTO names
        a[0] = a[0].replace(' ', '_')
        a[0] = a[0].replace('(', '')
        a[0] = a[0].replace(')', '')
        x = a[0]
        y = a[1]

        # some has and typed in for conflicts: 'O and 1' replace this with '0&1'
        if 'and' in y:
            y = y.replace(' and ','&')

        # store all the taxa in the input matrix in this list below
        alltaxa.append(x)

        # remaining species are stored in the list below

        if (y == '?'):
            missinglist.append(x)
        else:
            remainingspecies[x] =y

    # saving the character name/ column name in a variable
    if 'taxa_name' in line:
        line = line.strip()
        a = line.split('\t')
        char_name = a[1]



######################################################################################################################
#writing the missing removed matrix

out2 = open(os.path.join("intermediate_matrices","preprocessed_matrix.txt"), 'wb+')
# writing the correct header with the character name
out2.write('taxa_name\t%s\n'%(char_name))
for i in remainingspecies:
    out2.write('%s\t%s\n'%(i,remainingspecies[i]))

out2.close()

######################################################################################################################
# writing all the taxa in the original matrix in originaldatamatrix_taxalist file

out = open(os.path.join("statistics","originaldatamatrix_taxalist.txt"), 'wb+')

for i in alltaxa:
    out.write('%s\n'%(i))
print '\n'
print 'Step 2: Matrix pre-processing.......'
print 'removing the taxa with missing data from the matrix'
print 'Statistics for total taxa and missing character state taxa'
print 'Number of total taxa in the input matrix:',len(alltaxa)
print 'Number of taxa with missing data(?)  that was removed from the matrix',len(missinglist)
print 'Number of remaining taxa',len(remainingspecies)

out.close()
#print len(set(missinglist))

######################################################################################################################
######################################################################################################################
# writing the missing taxa separated into different levels

out1 = open(os.path.join("statistics","missingtaxa.txt"), 'wb+')

out1.write('total missing taxa: %s\n'%(len(missinglist)))

mt.taxaseparate(missinglist,out1,rank)

out1.close()

######################################################################################################################
# writing all taxa in the matrix separeted into levels
out3 = open(os.path.join("statistics","originaldatamatrix_taxalist_separated.txt"), 'wb+')
out3 = open('statistics/originaldatamatrix_taxalist_separated.txt', 'wb+')
out3.write('total  taxa: %s\n'%(len(alltaxa)))

mt.taxaseparate(alltaxa,out3,rank)

out3.close()