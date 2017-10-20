# this file works for the new format
#detects taxa with '0&1' and seperates them according to taxa level
# prints the literature sources for each taxa : the output is intented tobe reading only
# this code also replaces the intermediate taxa state of internal nodes by '?'
# 06/28/16


'''Instructions for running the code
This code must be run twice for each fin (pectoral and pelvic fin)
The reason is it requires the pectoral and pelvic xml data to identify the conflict literature sources
Make sure the pectoral fin xml file contains the pectoral name in the file name
You can save the intermediatecounts.txt file as pectoral and pelvic fin conflicting data'''


__author__ = 'pasan fernando'


import re

import collections
import sys
import os



sys.path.insert(0, os.path.join("tree_integration_pipeline","lib"))
import methods as mt


#printing the description of this code
print '\n'
print 'Step 3: Removing 0&1 states from higher-level taxa(families,genera,etc.)......... '



####################################################################################################################
# removing the internal taxa containing intermediate states
# note that when writing the new names the space in species name is replaced by '_" from herein


# a dictionary to store the taxa with '0&1' states
char_dic =[]

# defining the name of the output file
out1 = open(os.path.join("intermediate_matrices","conflicts_removed_datamatrix.txt"), 'wb+')

# writing the header in the output file


# Opening the input matrix and removing the conflict states '0&1's from only the interanal nodes
da = open(os.path.join("intermediate_matrices","preprocessed_matrix.txt"), 'r')

for line in da:
    #print line
    if (line != '\n'):
        line = line.strip()
        a = line.split('\t')
        a[0] = a[0].strip('\'')
        a[0] = a[0].replace(' ', '_')
        # saving the character name/ column name in a variable
        if 'taxa_name' in line:
            char_name = a[1]
            char_name = char_name.replace(' ', '_')

            # writing the header row in the output file. Since taxaname is only found once it will print the header only once
            out1.write('taxa_name\t%s\n' % (char_name))



        else:
            x = a[1]

            # the following if statement excludes all the species and removes '0&1' states only from internal nodes

            if a[1] == '0&1':
                #print a[0]
                # to check with the meta data xml file, the underscore needs to be replaced with a space
                char_dic.append(a[0].replace('_', ' '))
                if '_' not in a[0]:
                    x = '?'

            out1.write('%s\t%s\n'%(a[0],x))

out1.close()
print 'its running'
    ####################################################################################################################
    # if only the user has meta data file execute the following codes
def conflict_state_remove(metafile):
    # openinng the xml file from inputs folder


    p = open(os.path.join("inputs",metafile), 'r')


    # reading the names of each taxa and storing in a dictionary

    # This dictionary stores VTO id as key and taxa name as the value
    name = {}
    # a multiple dictionary for storing literature sources
    litsource = collections.defaultdict(list)
    # a counter for counting the number of taxa with '0&1' states
    #count = []
    # counter for counting the number of conflicting states
    count = 0

    '''the code section below, reads the xml file for the given fin, and checks wether they have '0&1' states in the matrix
       Then it searches for the literature sources in the xml file and writes the conflicting sources in intermediatecounts.txt file'''

    for line in p:

        if '<otu id=' in line:
            result = re.search('<otu id="(.*)" label', line)
            x = result.group(1)
            nm = re.search('label="(.*)" about', line)
            y = nm.group(1)
            name[x] = y


        if '<row id=' in line:
            result = re.search('otu="(.*)"', line)
            x = result.group(1)
            #print x
            # print name[x]

            if name[x] in char_dic:
                #print name[x]
                count = count +1


        # This code finds the literature source for each occurance of conflict state
        if '<meta xsi:type="LiteralMeta" property="dc:source">' in line:
            result = re.search('<meta xsi:type="LiteralMeta" property="dc:source">(.*)</meta>', line)
            y = result.group(1)
            #print y
            # store the literacture sources of conflicting states in a multiple dictionare
            # taxa name is the key and conflicting literature sources are values
            if y not in litsource[name[x]]:
                litsource[name[x]].append(y)




    out = open(os.path.join("statistics","conflict_counts.txt"), 'wb+')



    # writing conflict states for pectoral fin

    out.write('polymorpic state statistics for %s \n'%(char_name))
    out.write('total number of taxa: %s\n'%(len(char_dic)))
    out.write('\n')
    out.write('\n')
    mt.conflict_counter(char_dic,out,litsource)
    out.close()
    p.close()
    return



