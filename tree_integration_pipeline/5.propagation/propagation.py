
# Date 06/06/2016
# This scripts propagates data of internal nodes to their species.
# If the species already have data for a specific internal node that data will be kept.
# All the other species without data will be added to data file. For now, the propagation only considers internal taxa up to family level.
__author__ = 'pasan fernando'
######################################################################################################################

import re
import networkx as nx
import collections
import sys
import os

def propagation(meta_state):
    sys.path.insert(0, os.path.join("tree_integration_pipeline", "lib"))
    import methods as mt

    # printing the description of the code
    print '\n'
    print 'Step 5: Propagating data from families and genera....'
    ###### Reading the Vertebrate Taxonomy ontology (VTO) file to store relationships

    name,namer,rank,G = mt.VTOreader()


    # the dictionaries and lists required for propagation method is listed below
    character_dic = {}
    character_original = {}

    families =[]
    genus = []
    species = []
    order =[]
    nochild= []
    other =[]


    ######Read the data file separate the taxa into levels: family, order, species, genera, others (differentiate between genus and others in this step)
    if meta_state=='y':
        print 'you have a meta-data file; the output matrix will have 4 columns including the inferred state column'
        d = open(os.path.join("intermediate_matrices", "modified_inferredadded_matrix.txt"), 'r')

    else:
        print 'you do not have a meta-data file; the output matrix will have 3 columns excluding the inferred state column'
        d = open(os.path.join("intermediate_matrices", "conflicts_removed_datamatrix.txt"), 'r')

    # defining the propagated matrix; this is the final matrix that contains VTO names
    out = open(os.path.join("intermediate_matrices", "finalVTOmatrix.txt"), 'wb+')


    for line in d:
        if (line != '\n'):
            line = line.strip()
            a = line.split('\t')
            # saving the character name/ column name in a variable
            if 'taxa_name' in line:
                char_name = a[1]
                char_name = char_name.replace(' ', '_')

                # if statement is used to distinguish between user having meta data file or not. If user has the meta-data file, the output matrix contains
                # a column for inferred states, otherwise it does not.
                if meta_state=='y':

                    # writing the header row in the output file. Since taxaname is only found once it will print the header only once
                    out.write('taxa_name\t%s\t%s_inferred\t%s_propagated\n' % (char_name, char_name,char_name))

                else:
                    out.write('taxa_name\t%s\t%s_propagated\n' % (char_name, char_name))


            else:

                a[0] = a[0].strip('\'')
                # a[0] = a[0].replace(' ', '_')
                # a[0] = a[0].replace('(', '')
                # a[0] = a[0].replace(')', '')
                x = a[0]

                #if the user inputs the meta-data file, the propagation is performed on character states in the inferred column
                if meta_state == 'y':
                    y = a[2]

                # Else, propagation is performed on the original data column (no inferred data is available)
                else:
                    y = a[1]

                character_original[x]= a[1]# saving original character_dic in a dictionary

                character_dic[x] = y  #saving character data in the dictionaries (this can be inferred state, given the meta-data file)

                #separating the taxa into different taxonomic levels and storing in the lists defined above
                if '_' in x:
                    species.append(x)
                elif x.endswith('idae'):
                    families.append(x)
                elif x.endswith('iformes'):
                    order.append(x)
                else:
                    if x in rank:
                        if rank[x] == '0000005':
                            genus.append(x)
                        else:
                            other.append(x)
                    else:
                        nochild.append(x)  # taxonomic level undefined




    # the list below contains species list for original set of species in the input matrix for this code
    speciesold = species
    speciesold = set(speciesold)
    # print nochild
    # print genus
    # print other

    ############ propagation
    character_new =[]
    newspecies = []

    ## these 3 dics are used for count the number of species in interanal taxa
    excount ={}
    necount ={}
    totcount ={}


    #### The method for performing the propagation

    def propagate(x):
        replaced = []
        empty =[]
        character_replaced =[]



        # iterating trough internal taxa in the given list
        for i in x:

            # if taxa do not have data, append them to empty list
            if (character_dic[i]== '?'):
                empty.append(i)

            # if the taxa do have data...
            else:
                # append them to the replaced list
                replaced.append(i)

                # retrieving the descendents of each internal taxa
                li = list(nx.descendants(G, name[i]))
                tot=0
                ec =0
                nc =0

                # for each taxa in the descendents list
                for j in li:
                    # retrieve the taxa name
                    j = namer[j]

                    # check wether it is a species or genus
                    if ' ' in j:
                        tot=tot+1
                        j = j.replace(' ','_')

                        # if the species is in the input matrix go into this code
                        if j in species:
                            ec= ec+1
                            # propagate only if the existing species do not have data or have '?' as state
                            # propagating the character
                            if character_dic[j]== '?':
                                if character_dic[i] != '?':
                                    character_dic[j] = character_dic[i]
                                    character_new.append(j)
                                    character_replaced.append(i)


                        # if the descendent species is not in the original matrix we need to add it
                        else:
                            # counting the newly added species
                            nc = nc+1
                            # appending new species into a list
                            newspecies.append(j)
                            species.append(j)
                            # propagating the data into the character
                            character_dic[j] = character_dic[i]

                            # counting the propagated number for newly added species for the character
                            if character_dic[i] != '?':
                                character_new.append(j)
                                character_replaced.append(i)

                excount[i]=ec
                necount[i]=nc
                totcount[i]=tot
        # the method peforms the propagation and retuns lists for empty taxa, total propagated, and character propagated
        return replaced,empty,character_replaced


    # these two lines actually performs the propagation

    # first start propagating the genus data; this must be peformed before family data propagation
    genus_replaced, genus_empty, genus_pec = propagate(genus)

    # then propagating the family data
    family_replaced, family_empty, family_pec = propagate(families)

    ########### printing the final matrix output




    # The propagated matrix only contains species; iterate through species

    for i in species:
        #print i
        # here, pc and pl are the original character states for each species; for instance, state 2 is inferred presence
        # which is originally represented as 1

        if character_dic[i] =='2':
            pc = '1'

        # If the inferred state is '3' (inferred absence), the original state is '0'
        elif character_dic[i] =='3':
            pc = '0'

        # else, it is asserted: 0 or 1
        else:
            pc = character_dic[i]


        # if the species was propagated for the character fin, it should be found within character_new list, thus represented as 1 for propagation
        if i in character_new:
            pcn = '1'
        else:
            pcn = '0'

        # Again, if meta-data is available the output matrix contains four columns, including the inferred state column
        if meta_state == 'y':
            out.write('%s\t%s\t%s\t%s\n'%(i,pc,character_dic[i],pcn))

        # else, the output matrix contains only 3 columns, excluding the inferred state column
        else:
            out.write('%s\t%s\t%s\n' % (i, pc, pcn))

    out.close()

    ###################### printing the statistics#######################################

    out1 = open(os.path.join("statistics", "propagationstatistics.txt"), 'wb+')

    out1.write('************* statistics for the original data file *******\n')
    out1.write('number of species: %s\n'%(len(speciesold)))
    out1.write('number of genera: %s\n'%(len(genus)))
    out1.write('number of families: %s\n'%(len(families)))
    out1.write('number of orders: %s\n'%(len(order)))
    out1.write('number of other higher level taxa: %s\n'%(len(nochild)))


    out1.write('\n')
    out1.write('\n')





    out1.write('************* statistics for the new data file *******\n')
    print 'propagation is complete!'
    print len(newspecies),'new species were added to the matrix due to propagation'
    out1.write('number of newly added species: %s\n'%(len(newspecies)))
    print len(species),'species are total in the final matrix'
    out1.write('number of total species in the new data file: %s\n'%(len(species)))

    out1.write('\n')
    out1.write('\n')

    out1.write('************* propagation statistics *******\n')

    out1.write('newly propagated species for %s: %s\n'%(char_name,len(set(character_new))))
    genus_replaced = set(genus_replaced)
    family_replaced = set(family_replaced)
    out1.write('number of genera with data: %s\n'%(len(genus_replaced)))
    out1.write('number of families with data: %s\n'%(len(family_replaced)))
    out1.write('\n')
    genus_pec = set(genus_pec)
    family_pec = set(family_pec)
    out1.write('Number of families that propagated data for %s: %s\n'%(char_name,len(family_pec)))
    out1.write('Number of genera that propagated data for %s: %s\n'%(char_name,len(genus_pec)))
    out1.write('\n')
    out1.write('\n')



    out1.write('there were some internal nodes that had ? for the character. There is no use of them for the propagation. They are printed below\n')
    out1.write('\n')
    out1.write('families :')
    for i in family_empty:
        out1.write('%s,'%(i))
    out1.write('\n')
    out1.write('\n')

    out1.write('genera :')
    for i in genus_empty:
        out1.write('%s,'%(i))

    out1.write('\n')
    out1.write('\n')

    out1.write('******the propagated counts for each internal node is given below\n')
    out1.write('\n')
    out1.write('families\t#existing_species_count\tnew_species_count\ttotal_number_of_speciesinVTO\n')
    famcount=0
    for i in family_replaced:
        famcount=famcount+totcount[i]
        out1.write('%s\t%s\t%s\t%s\n'%(i,excount[i],necount[i],totcount[i]))
    out1.write('\n')
    out1.write('genus\t#existing_species_count\tnew_species_count\ttotal_number_of_speciesinVTO\n')
    for i in genus_replaced:
        out1.write('%s\t%s\t%s\t%s\n'%(i,excount[i],necount[i],totcount[i]))

    out1.close()
    ######printing a list of all species taxa for comparison purposes

    out2 = open(os.path.join("statistics", "finalVTOspecieslist.txt"), 'wb+')

    for i in species:
        out2.write('%s\n'%(i))

    out2.close()

    #print excount
    #print necount
    #print famcount
    #print sum(totcount.values())
    #print sum(necount.values())



    # end of the propagation method
    return