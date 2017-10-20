import dendropy
import re
import os

print "Please move the data matrix (nexml file) downloaded from Phenoscape KB to the inputs folder"



#opening the xml file
while True:
    filename = raw_input('Please enter the name of the nexml file downloaded from phenoscape:').lower()

    # nexmlname = 'inputs/'+filename
    nexmlname = os.path.join("inputs", filename)

    # nexmlname = "../inputs/teleostei-pectoralfin.xml"

    # Using error handling to prompt the user to correctly enter the file name
    try:
        xmlfile = open(nexmlname, 'r')
        break
    except (IOError):
        print "The file name you entered is wrong. It is not found in the inputs folder. Please retry!"
        continue



# defining empty list to store character names
charlist = []

#this loop appends character names to the empty list
for line in xmlfile:
    if '<char id=' in line:
        #print line
        name1 = re.search('label="(.*)" about=', line)
        cname = name1.group(1)
        #print cname
        charlist.append(cname)

# printing the character list
print '\n'
print 'Step1: converting the nexml matrix to tab-delimited format........'
print 'The character in your matrix is:',cname

# converting the character list into a tuple
chartuple = tuple(charlist)

#print chartuple

# reading the nexml matrix using dendropy module
cmatrix = dendropy.StandardCharacterMatrix.get(path=nexmlname,schema="nexml")

# opening the output matrix
out = open("intermediate_matrices/tabdelemited_charactermatrix.txt", 'wb+')

#writing the header of the output matrix: first the taxa name
out.write('taxa_name\t')

# and then all the character names
for i in chartuple:
    out.write('%s\t'%(i))

out.write('\n')

tot_num_char = len(chartuple)

# following loop goes fills character states for all taxa in tabdelimited format
# for conflicts it will print as 0 and 1, which is the default representation in phenoscape
# all the empty cells (without data) will have ?
for taxon in cmatrix:
    charlist = cmatrix[taxon]
    out.write('%s\t'%(taxon.label))
    #out.write('%s\n'%(charlist))
    number_to_append = tot_num_char - len(charlist)
    for num in range(number_to_append):
        charlist.append(None)
    for i in charlist:
        if i == None:
            i = '?'
        out.write('%s\t' % (i))
    out.write('\n')

out.close()