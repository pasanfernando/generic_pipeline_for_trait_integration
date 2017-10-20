
# Date 06/12/2016
#
# This script is used to convert the VTO data matrix into open tree version.
# It will use combined approach that includes name conversion using NCBI ids and direct name matching.
# The script can also generate statistics for the final conversion.
__author__ = 'pasan fernando'
######################################################################################################################


import sys
import re
import networkx as nx
import collections
import dendropy
import os

sys.path.insert(0, os.path.join("tree_integration_pipeline", "lib"))
import methods as mt





def taxonomic_reconciliation(metastate):

    #printing the code description
    print '\n'
    print 'Step 6: this is the last step of the pipeline; The algorithm will try to match taxa names of Phenoscape data matrix'
    print 'with the taxa names of open tree using both NCBI taxon ids and taxonomic names'




    ############ Reading the Vertebrate Taxonomy ontology (VTO) file to store relationships##############################################################################
    # calling the advanced VTO reader method from the methods library



    name,namer,vtncbi,vtncbir,extinct,G =mt.advancedVTOreader()

    ########################## open tree data file processing########################################################################

    otncbi ={}
    otncbir={}
    oid ={}

    # opening the taxonomy file downloaded from opentree; this contains references to NCBI ids
    d = open(os.path.join("inputs", "taxonomy.tsv"), 'r')


    for line in d:
        a = line.split('|')

        # extracting open tree id and name
        otid = a[0].strip()
        otname = a[2].strip()

        #storing all ids: name as the key and open tree id as the value
        oid[otname]=otid

        # extracting the NCBI ids, by stripping unwanted segments
        a1 = a[4].strip()

        if 'ncbi:' in a1:
            if ',' in a1:
                b = a1.split(',')
                for i in b:
                    if 'ncbi' in i:
                        x = i

            else:
                x = a1
            #print '     ',x
            x=x.strip('ncbi:')

            # storing ncbi ids: name as key ncbi ID as value
            otncbi[otname]=x
            #reversed assignment: key is ncbi id , value is name
            otncbir[x]=otname

    ##### processing the VTO matrix############################################################################################################
    matchnames={}
    matchnamesr ={}
    vtlist=[]
    vtwitnc =[]
    vtdic ={}

    # list to count absences
    peclosslist =[]
    pecprop =[]
    pec ={}

    char_name= ''

    # read the final VTO matrix, separate it so that name is the key and the rest of the line is the value
    da = open(os.path.join("intermediate_matrices", "finalVTOmatrix.txt"), 'r')

    for line in da:
        line = line.strip()
        b = line.split('\t')
        a = line.partition('\t')
        # the if statement here is used to extract the species with pectoral or pelvic fin absences
        # we are interested in knowing how many absencses in the final matrix are mismatched and propagated
        if (b[1] =='0') or (b[1] =='0&1'):
            # replacing the underscores of taxa names by space

            a1 = a[0].replace('_', ' ')
            a1 = a1.strip('\"')


            peclosslist.append(a1)
            pec[a1]=b[1]
            # extracting the propagated absences for pectoral fin
            if metastate == 'y':
                if b[3] =='1':
                    pecprop.append(a1)

            else:
                if b[2] =='1':
                    pecprop.append(a1)




        # this code works for all the taxa; not only the absences
        if a[0] =='taxa_name':
            char_name = b[1]
            char_name = char_name.replace(' ', '_')

        else:
            a1= a[0].replace('_',' ')
            a1= a1.strip('\"')
            vtlist.append(a1)
            vtdic[a1]=a[2]


    print 'the number of propagated VTO data matrix taxa:',len(vtlist)

    #print char_name
    ####################################################################################################################################
    # defining a list to store VTO names with ncbi ids that is matching in OT
    vtwitncid =[]
    # defining the list to store VTO taxa that match by name and another list for name mismatches
    otmatched =[]
    otmismatched =[]


    ## checking the VTO list taxa has ncbi id and in Open Tree

    for i in vtlist:

        if i in vtncbi:

            nc1 =vtncbi[i]
            if nc1 in otncbir:
                vtwitncid.append(i)
                nm1 =otncbir[nc1]
                matchnames[i]=nm1
                #reversed dic : key is open tree name
                matchnamesr[nm1]=i
        else:
            vtwitnc.append(i)

    #if it does not have ncbi id append to another list and check wether they have the name in open tree
    for i in vtwitnc:
        if i in oid:
            matchnames[i]=i
            matchnamesr[i]=i


    ################################################################################################################################################
    # reading the open tree taxonomy file in newick format and storing the taxa names in a list
    otlist=[]
    # ot = open('all_tips.txt', 'r')
    # for line in ot:
    #     a1 = line.strip()
    #     otlist.append(a1)
    print 'This step requires the open tree file you downloaded from Open Tree of Life'
    while True:
        treefile = raw_input('Please move the open tree file to input folder and type its name with the extension:').lower()

        # Using error handling to prompt the user to correctly enter the file name
        try:
            treepath = os.path.join("inputs", treefile)
            tree1 = dendropy.Tree.get(path=treepath, schema="newick", preserve_underscores=True)
            break

        except (IOError):
            print "The file name you entered is wrong. It is not found in the inputs folder. Please retry!"
            continue


    #print tree1
    for i in tree1:
        #print i.label
        #print i.taxon
        tname= str(i.taxon)
        tname=tname.strip('\'')
        if tname != 'None':
            # this code removes the ott header from the end of the
            tname = re.sub(r'(_| )ott\d{3,8}', '', tname)


            # replacing the underscores in the name by a space
            tname = tname.replace('_',' ')
            otlist.append(tname)
            # writing a ott removed newick file so the data matrix matches with tree file

            taxa1 = dendropy.datamodel.taxonmodel.Taxon(label=tname)
            i.taxon = taxa1

    # writing a ott removed newick file so the data matrix matches with tree file
    tree1.write(path="outputs/preprocessed_tree.tre", schema="newick", unquoted_underscores=True)

    # writing the final OT taxa list in  a file
    out3 = open(os.path.join("statistics", "opentree_specieslist.txt"), 'wb+')

    for i in otlist:
        out3.write('%s\n'%(i))
    out3.close()
    ################################################################################################################################################

    ### matching with the final tree list
    # checking the species that is not matching by name but that has ncbi ids
    for i in vtlist:
        if i in otlist:
            otmatched.append(i)
        else:
            otmismatched.append(i)

    matched = matchnames.keys()

    # extracting the list of mismatched species
    mismatched = set(vtlist) - set(matched)


    matchedvals = matchnames.values()

    # print 'VTO name count without ncbi ID:', len(vtwitnc)
    # print 'mapped count from VTO to OT conversion (using database):', len(matchnames)
    # print 'final mismatched list (database comparison): ', len(mismatched)

    # Here, we use the database comparison. The reason is that some times some of the taxa in VTO matrix cannot be found in
    # Open Tree phylogeny, but can be found in Open Tree data base.

    # getting the initial database vs tree file mismatched
    fmismatch = set(matchedvals) - set(otlist)

    # print fmismatch

    ## the below code was required because there are name changes between OT database file and OT tree file
    # during the initial matching
    # checkout Danio rerio for more information in the OT data file and the tips list

    ## so this is the second round of matched dic replacement
    #print 'initial mismatch count between database and the final tree file', len(fmismatch)
    for i in fmismatch:
        k = matchnamesr[i]

        if k in otlist:
            del matchnamesr[i]
            # print k
            matchnames[k] = k
            matchnamesr[k] = k

            # print mismatched
            # print vtwitnc
    # the intersection of name mismatches with ncbi id has ones
    contranames = set(otmismatched) & set(vtwitncid)

    ### generating second round statistics##############
    matchedot1 =matchnamesr.keys()
    ffmismatch = set(matchedot1)- set(otlist)
    # print 'Final mismatch count between database and the final tree file',len(ffmismatch)
    # print len(fmismatch)-len(ffmismatch), 'mismatches were solved by second round'



    #################### generating the final open tree matrix####################################################################################

    # defining the output file for full open tree species. Including the ones without data
    out1 = open(os.path.join("outputs", "finalfullopentree_matrix.txt"), 'wb+')


    # defining the output file for only the matched species with data
    out2 = open(os.path.join("outputs", "finalopentree_matrix_onlydata.txt"), 'wb+')


    # writing the header in each output file (using char_name variable)
    if metastate == 'y': # if meta-data is available, 4 columns are in the output matrix including inferred column
        out1.write('taxa_name\t%s\t%s_inferred\t%s_propagated\n' % (char_name, char_name,char_name))
        out2.write('taxa_name\t%s\t%s_inferred\t%s_propagated\n' % (char_name, char_name,char_name))
    else: # else, 3 columns are in the output matrix excluding inferred column
        out1.write('taxa_name\t%s\t%s_propagated\n' % (char_name, char_name))
        out2.write('taxa_name\t%s\t%s_propagated\n' % (char_name, char_name))

    # read the opentree tips list if the taxa is in the mapped dictionary
    #   get the data from the vtodic and print it out, else: print question marks
    matchedots =[]
    matchedvts =[]

    for i in otlist:
        if i in matchnamesr:
            matchedots.append(i)
            a = matchnamesr[i]
            matchedvts.append(a)
            out1.write('%s\t%s\n'%(i,vtdic[a]))
            out2.write('%s\t%s\n' % (i, vtdic[a]))
        else:
            if metastate == 'y': # if meta-data is available, 4 columns are in the output matrix including inferred column
                out1.write('%s\t?\t?\t?\n'%(i))
            else:# else, 3 columns are in the output matrix excluding inferred column
                out1.write('%s\t?\t?\n' % (i))

    out1.close()
    out2.close()


    # printing the statistics
    # print 'matched count of OT ids: ', len(matchedots)
    # print 'matched count of VTO ids: ', len(matchedvts)
    mismatchedvts = set(vtlist) - set(matchedvts)
    mismatchedvts=list(mismatchedvts)
    # print 'final mismatched count in VTO', len(vtlist)-len(matchedvts)
    # print 'final mismatched count in VTO', len(mismatchedvts)

    # seperating the mismatched VTO list

    spmismatches =[]

    #Detecting extinct taxa

    extinctmismatches = set(mismatchedvts)& set(extinct)
    ermismatchedvts = set(mismatchedvts) - extinctmismatches

    #detecting sp.
    #
    for i in ermismatchedvts:
        if ('sp.'in i) or ('cf.'in i) or ('Species'in i) or ('Genus'in i):
            spmismatches.append(i)

    finalmismatchedvts = ermismatchedvts - set(spmismatches)

    ################################################################################################################################################
    #update: 11/2/2016: this snippet was added to print the taxa that do not have ncbi ids and the ones that have them and do not match with OT
    # total set of mismatches based only on ncbi id matching
    ncidmismatches= set(vtlist)-set(vtwitncid)

    #getting the ones that have NCBI id but do not match with OT

    otncidmismatch=ncidmismatches-set(vtwitnc)
    #print 'ones that have NCBI id but do not match with OT', len(otncidmismatch)


    ################################################################################################################################################
    ## printing the mismatched list
    out = open(os.path.join("statistics", "finalmismatchedlist_andstats.txt"), 'wb+')

    out.write('VTO data matrix taxa count:%s\n'%(len(vtlist)))
    out.write('Open tree file taxa count:%s\n'%(len(otlist)))
    print 'Number of taxa in the open tree phylogeny: ', len(otlist)

    out.write('\n')
    out.write('The number of VTO taxa that is mismatched with OT tree file:%s\n'%(len(mismatchedvts)))

    out.write('The number of VTO taxa that is matched with OT tree file:%s\n'%(len(matchedvts)))
    print 'This algorithm managed to match',len(matchedvts),'taxa with the open tree file'
    out.write('\n')

    out.write('####################################################################################\n')

    out.write('####################################################################################\n')
    out.write('The number of mismatched taxa that is extinct:%s\n'%(len(extinctmismatches)))
    for i in extinctmismatches:
        out.write('%s\n'%(i))

    out.write('\n')
    out.write('\n')

    out.write('The number of mismatched taxa that has sp. in name (improper naming):%s\n'%(len(spmismatches)))
    for i in spmismatches:
        out.write('%s\n'%(i))

    out.write('\n')
    out.write('\n')

    out.write('The remaining mismatches list:%s\n'%(len(finalmismatchedvts)))
    for i in finalmismatchedvts:
        out.write('%s\n'%(i))

    out.write('\n')
    out.write('\n')

    out.write('out of the remaining mismatched taxa There are %s VTOtaxa that is in the Open tree Database but they were not in the final open tree file\n'%(len(ffmismatch)))
    out.write('VTO_name\tOTname\n')
    for i in ffmismatch:

        out.write('%s\t%s\n'%(matchnamesr[i],i))
    out.write('\n')
    out.write('\n')
    # out.write('There are %s VTOtaxa that is not matched in Open tree Database\n'%(len(mismatched)))
    # out.write('\n')
    # for i in mismatched:
    #     out.write('%s\n'%(i))
    out.write('####################################################################################\n')
    out.write('the mismatched taxa with %s loss\n'%(char_name))
    out.write('\n')

    pecpropmismatchloss = set(peclosslist)&set(mismatchedvts)&set(pecprop)

    pecunpropmismatchloss = (set(peclosslist)&set(mismatchedvts))-set(pecprop)


    out.write('the mismatched taxa for %s propagated with loss: %s\n'%(char_name,len(pecpropmismatchloss)))
    for i in pecpropmismatchloss:
        out.write('%s\t%s\n' % (i,pec[i]))
    out.write('\n')

    out.write('the mismatched taxa for %s unpropagated with loss: %s\n'%(char_name,len(pecunpropmismatchloss)))
    for i in pecunpropmismatchloss:
        out.write('%s\t%s\n' % (i,pec[i]))
    out.write('\n')



    out.write('####################################################################################\n')
    out.write('There are %s that that match by only NCBI ids \n'%(len(vtwitncid)))
    out.write('There are %s that that is mismatched by only NCBI ids \n'%(len(vtlist)-len(vtwitncid)))
    out.write('\n')
    out.write('\n')
    out.write('There are %s that that match by only names \n'%(len(otmatched)))
    out.write('There are %s that that is mismatched by only names \n'%(len(otmismatched)))
    out.write('There are %s that do not match by name but are matching by NCBI ids \n'%(len(contranames)))
    out.write('\n')
    out.write('VTO_name\tOTname\n')
    for i in contranames:
        out.write('%s\t%s\n' % (i,matchnames[i]))

    out.write('\n')
    out.write('There are %s that that do not have ncbi ids in vt \n'%(len(vtwitnc)))
    for i in vtwitnc:
        out.write('%s\n'%(i))
    out.write('\n')
    out.write('There are %s that that do  have ncbi ids in vt but do not match with OT \n' % (len(otncidmismatch)))
    for i in otncidmismatch:
        out.write('%s\n'%(i))

    out.close()

    ################################################################################################################################################
    # writing the final list of matched taxa from data matrix
    out4 = open(os.path.join("statistics", "finalmatchedvtlist.txt"), 'wb+')


    for i in matchedvts:
        out4.write('%s\n'%(i))
    out4.close()

    return


