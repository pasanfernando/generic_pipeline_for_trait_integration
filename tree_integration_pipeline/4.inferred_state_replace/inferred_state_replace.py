# author: Pasan Fernando
# Date: 06/15/16
# Used to replace the character state of data matrix file for the taxa with inferred presence and absence

#################################################################################################

# defining lists to store character inferred presence
pec = []
import re
import collections
import sys
import os

''' To run this code, taxa with inferred presence and absence must be extracted from meta-data. Hence, the metadata file entered when running the main script will be 
passed to method in this script'''

print '\n'
print 'Step 4:Distinguishing inferred vs asserted states based on the information from meta-data........'

def inferred_state_replacer(metafile):

    sys.path.insert(0, os.path.join("tree_integration_pipeline", "lib"))
    import methods as mt

    # openinng the xml file from inputs folder
    xmlfile = open(os.path.join("inputs", metafile), 'r')


    # defining empty list to store character names
    charlist = []
    #reading the names of each taxa
    name = {}


    ontab = collections.defaultdict(list)
    ontpr = collections.defaultdict(list)
    count = []

    s = None
    #this loop appends character names to the empty list
    for line in xmlfile:
        #print line
        if '<char id=' in line:
            #print line
            name1 = re.search('<char id="(.*)" label=', line)
            cid = name1.group(1)
            #print cname
            charlist.append(cid)

            if len(charlist) == 1:
                charid = charlist[0]
                #print charid
                s = '<meta xsi:type="LiteralMeta" property="dc:identifier">' + charid
            else:
                print "There are more than one characters in your metadata file. There should be only one character\n"

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

        if s != None:
            if s+'_0' in line:

                z = '0'

            if s+'_1' in line:
                z = '1'

        if '<meta xsi:type="LiteralMeta" property="ps:isDirect">' in line:
            #print 'yes'
            result = re.search('<meta xsi:type="LiteralMeta" property="ps:isDirect">(.*)</meta>', line)
            y = result.group(1)
            # print y
            if z == '0':
                ontab[x].append(y)
            if z == '1':
                ontpr[x].append(y)

            if y == 'false':
                if z =='0':
                    count.append(x)

    #print ontab
    #print ontpr

    # the lists to sepearate asserted vs inferred presence vs absence

    assertedab = []
    inferredab = []
    assertedpr = []
    inferredpr = []

    assertedab,inferredab = mt.infervsasserted(ontab,name)
    assertedpr,inferredpr = mt.infervsasserted(ontpr,name)

    #print 'inferred absence:',inferredab
    #lists to store matched inferred presence and absence taxa
    matchedinfpr=[]
    matchedinfab=[]

    # reading the input data matrix for the code
    m = open(os.path.join("intermediate_matrices", "conflicts_removed_datamatrix.txt"), 'r')


    # defining the output matrix of the code
    out = open(os.path.join("intermediate_matrices", "modified_inferredadded_matrix.txt"), 'wb+')

    # reading the input and selects the states with inferred presence which are represented by '2' hereafter
    # also if there any taxa with inferred absence only they will be represented by '3' hereafter

    for line in m:
        if (line != '\n'):
            line = line.strip()
            a = line.split('\t')
            a[0] =a[0].strip('\'')
            a1 = a[0]  # this line was required to keep the parenthesis within taxa names in the final VTO matrix
            a[0] = a[0].replace(' ', '_')

            # saving the character name/ column name in a variable
            if 'taxa_name' in line:

                char_name = a[1]
                char_name = char_name.replace(' ', '_')

                # writing the header row in the output file. Since taxaname is only found once it will print the header only once
                out.write('taxa_name\t%s\t%s_inferred\n' % (char_name,char_name))

            else:
                a[0] = a[0].replace('(', '')
                a[0] = a[0].replace(')', '')
                if a[0] in inferredpr:
                    # checking wether there are no asserted absences; there can be when there are '0&1' states
                    if a[0] not in assertedab:
                        x = '2'
                        matchedinfpr.append(a[0])
                    # if there are asserted absences keep the original '0&1' state
                    else:
                        x = a[1]

                elif a[0] in inferredab:
                    # checking wether there are no asserted presence; there can be when there are '0&1' states
                    if a[0] not in assertedpr:
                        x = '3'
                        matchedinfab.append(a[0])
                    # if there are asserted presences, keep the original '0&1' state
                    else:
                        x = a[1]
                else:
                    x = a[1]



                out.write('%s\t%s\t%s\n'%(a1,a[1],x))
    out.close()
    #################################################################################################
    # generating statistics for inferred data

    # defining another output file for inferred state statistics
    out1 = open(os.path.join("statistics", "inferredstats.txt"), 'wb+')


    # writing the number of taxa with inferred presence state
    #out1.write('inferred presence for pectoral fin: %i\n'%len(inferredpr))

    #writing a list of taxa seperated into differenet levels

    out1.write('for asserted absence data\n')
    mt.inferred_taxa_separate(assertedab,out1)
    # pasan(empty)
    out1.write('\n')
    out1.write('for inferred absence data\n')
    mt.inferred_taxa_separate(matchedinfab,out1)
    out1.write('\n')

    out1.write('for asserted presence data\n')
    mt.inferred_taxa_separate(assertedpr,out1)
    # pasan(empty)
    out1.write('\n')
    out1.write('for inferred precence data\n')
    mt.inferred_taxa_separate(matchedinfpr,out1)
    out1.write('\n')

    # Counting the unmapped taxa with inferred presence; the unmapped taxa are mismatched between pectoral and pelvic xml files
    #with the input matrix with various naming errors; usually there are only one or two taxa like this; they are due to the presence
    # of quot instead of actual "" in the taxa name



    unmappedpr = set(inferredpr) - set(matchedinfpr)
    unmappedpr = list(unmappedpr)

    out1.write('unmapped data for pectoral fin presence: %i\n'%len(unmappedpr))
    for i in unmappedpr:
        out1.write('%s\n'%i)

    # repeating the same for unmapped absence

    unmappedab = set(inferredab) - set(matchedinfab)
    unmappedab = list(unmappedab)

    out1.write('unmapped data for pectoral fin absence: %i\n'%len(unmappedab))
    for i in unmappedab:
        out1.write('%s\n'%i)


    #print inferredpr

    out1.close()

    # printing some statistics on the screen
    print len(matchedinfpr),'inferred presence states and', len(matchedinfab), 'inferred absence states were replaced during this step'

    return

