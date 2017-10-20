
# Date 06/06/2016
# This scripts propagates data of internal nodes to their species.
# If the species already have data for a specific internal node that data will be kept.
# All the other species without data will be added to data file. For now, the propagation only considers internal taxa up to family level.
__author__ = 'pasan fernando'
######################################################################################################################

import re
import networkx as nx
import collections
import os

###### Reading the Vertebrate Taxonomy ontology (VTO) file to store relationships
#VTO reader method

def VTOreader():
    p = open(os.path.join("inputs", "vtonewfinal.owl"), 'r')


    G = nx.DiGraph()
    # The dictiorary to store VTO name as the key and id as the value
    name = {}
    # The dictiorary to store id as the key and VTO name as the value
    namer = {}
    # The dictionary to store the rank of the taxa name; rank changes based on taxnonomic level
    rank = {}

    # The following loop iterates the VTO.owl file and stores taxonomy hirachy in dictionaries
    for line in p:
        # print line
        if '<!-- http://purl.obolibrary.org/obo/' in line:
            result = re.search('<!-- http://purl.obolibrary.org/obo/(.*)-->', line)
            x = result.group(1)
            x1 = x.strip()
            # print x1
            if G.has_node(x1) == False:
                G.add_node(x1)

        if '<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">' in line:
            name1 = re.search('<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">(.*)</rdfs:label>', line)
            n = name1.group(1)
            # print n
            name[n] = x1
            namer[x1] = n
        #
        if '<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/' in line:
            s = re.search('<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/(.*)"/>', line)
            k = s.group(1)

            G.add_edge(k, x1)

        if '<vto:has_rank rdf:resource="http://purl.obolibrary.org/obo/TAXRANK_' in line:
            s = re.search('<vto:has_rank rdf:resource="http://purl.obolibrary.org/obo/TAXRANK_(.*)"/>', line)
            r = s.group(1)
            rank[n]=r

    return name,namer,rank,G

######################################################################################################################
# this is the advanced method for reading VTO taxa data
# used in the taxonomic_reconcilliation step
# in addition to the simple method, this makes ncbi to vto mapping dictionary and a list to store extinct taxa
def advancedVTOreader():

    p = open(os.path.join("inputs", "vtonewfinal.owl"), 'r')

    G = nx.DiGraph()
    # The dictiorary to store VTO name as the key and id as the value
    name = {}
    # The dictiorary to store VTO id as the key and VTO name as the value
    namer = {}
    # The dictiorary to store VTO id as the key and NCBI id as the value
    vtncbi = {}
    # The dictiorary to store NCBI id  as the key and VTO id as the value
    vtncbir = {}
    # a list to store extinct taxa
    extinct = []

    # The following loop iterates the VTO.owl file and stores taxonomy hierarchy in dictionaries
    for line in p:
        # print line
        if '<!-- http://purl.obolibrary.org/obo/' in line:
            result = re.search('<!-- http://purl.obolibrary.org/obo/(.*)-->', line)
            x = result.group(1)
            x1 = x.strip()
            # print x1
            if G.has_node(x1) == False:
                G.add_node(x1)

        if '<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">' in line:
            name1 = re.search('<rdfs:label rdf:datatype="http://www.w3.org/2001/XMLSchema#string">(.*)</rdfs:label>',
                              line)
            n = name1.group(1)
            # print n
            name[n] = x1
            namer[x1] = n
        #
        if '<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/' in line:
            s = re.search('<rdfs:subClassOf rdf:resource="http://purl.obolibrary.org/obo/(.*)"/>', line)
            k = s.group(1)

            G.add_edge(k, x1)
            # extracting the ncbi taxon
        if '<oboInOwl:hasDbXref rdf:datatype="http://www.w3.org/2001/XMLSchema#string">NCBITaxon:' in line:
            s = re.search('<oboInOwl:hasDbXref rdf:datatype="http://www.w3.org/2001/XMLSchema#string">NCBITaxon:(.*)<',
                          line)
            nc = s.group(1)
            ## there are synonyms for some taxa and they have ncbi id's as well so only the first ocurrance is stored
            if n not in vtncbi:
                vtncbi[n] = nc
            # reversed version: key is ncbi id
            if nc not in vtncbir:
                vtncbir[nc] = n

        # detecting and storing extinct taxa in VTO
        if '<vto:is_extinct rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">true</vto:is_extinct>' in line:
            extinct.append(n)
    return name,namer,vtncbi,vtncbir,extinct,G


######################################################################################################################
## a method to print all taxa in any given list in a given output file
def printout(lis, out):
    out.write('taxa count: %s\n'%(len(lis)))

    for l in lis:
        out.write('%s\n'%(l))

    return
######################################################################################################################
# a method to separate a list of taxa into different taxanomic levels basaed on VTO taxonomy ontology
# the separated taxa is written in another output file which is represented by the out parameter
def taxaseparate(lis1,out,rank):
    total = []
    families = []
    genus = []
    species = []
    order = []
    nochild = []
    other = []

    for x in lis1:
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
                nochild.append(x)

    out.write('orders\t')
    printout(order,out)
    out.write('\n')
    out.write('familes\t')
    printout(families,out)
    out.write('\n')
    out.write('\n')
    out.write('genera\t')
    printout(genus,out)
    out.write('\n')
    out.write('species\t')
    printout(species, out)
    out.write('\n')
    out.write('higher level taxa with rank\t')
    printout(other,out)
    out.write('\n')
    out.write('higher level taxa without rank\t')
    printout(nochild,out)
    out.write('\n')

    return

######################################################################################################################
# defining method for family count
# this method will be used by conflict_state_remover.py code
# This method takes a taxa list as an input and and separates them based on taxonomic level and writes the conflicting sources
# in output file
def conflict_counter(li,out,source):
    family = []
    genus = []
    species = []
    for e in li:
        nm1 = e
        if 'idae' in nm1:
            family.append(nm1)
        elif ' ' in nm1:
            species.append(nm1)
        else:
            genus.append(nm1)
    # print family
    # print species
    # print genus
    out.write('number of families: %d\n' % len(family))
    for line in family:
        li = source[line]
        out.write('%s\n' % (line,))
        if li:
            for i in li:
                out.write('       %s\n' % (i))
                # else:
                #     out.write('%s \n' % (line))
    out.write('\n')
    out.write('number of genera: %d\n' % len(genus))
    for line in genus:
        li = source[line]
        out.write('%s\n' % (line,))
        if li:
            for i in li:
                out.write('       %s\n' % (i))
    out.write('\n')
    out.write('number of species: %d\n' % len(species))
    for line in species:
        li = source[line]
        out.write('%s\n' % (line,))
        if li:
            for i in li:
                out.write('       %s\n' % (i))
    return


######################################################################################################################
# defining method for family count
# this method will be used by inferred_state_replace.py code
# This method takes a taxa list as an input and and separates them based on taxonomic level
# in output file


def inferred_taxa_separate(li,out):
    family =[]
    genus = []
    species = []
    for e in li:
        nm1 = e
        if '_' in nm1:
            species.append(nm1)
        elif 'idae' in nm1:
            family.append(nm1)

        else:
            genus.append(nm1)
    # print family
    # print species
    # print genus
    out.write('number of families: %d\n'%len(family))
    for line in family:
        out.write('%s\n'%(line))
    out.write('\n')
    out.write('number of genera: %d\n'%len(genus))
    for line in genus:
        out.write('%s\n'%(line))
    out.write('\n')
    out.write('number of species: %d\n'%len(species))
    for line in species:
        out.write('%s\n'%(line))
    return

######################################################################################################################
# this method will be used by inferred_state_replace.py code
# separating the absences or presences into asserted vs inferred
# If there is asserted evidence we put them into asserted list; do not care if they have inferred evidence
# only selected as inferred when there is no asserted evidence
# also the VTO ids are converted to names using the name dictionary and then the names are pre-processed accordingly

def infervsasserted(taxalist,name):
    asserted =[]
    inferred =[]
    for l in taxalist:
        x3 = taxalist[l]
        line = name[l]
        line = line.strip()
        line = line.replace(' ', '_')
        line = line.replace('(', '')
        line = line.replace(')', '')

        # print x3
        if 'true' in x3:
            asserted.append(line)
        else:
            inferred.append(line)

    return asserted,inferred
