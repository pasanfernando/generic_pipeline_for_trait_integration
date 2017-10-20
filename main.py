
import sys
import os

'''This is the main script that runs the pipeline, which converts the data matrix (supermatrix) downloaded from Phenoscape to a version that can be
merged with an Open Tree phylogeny. The pipeline has 6 basic steps. Some of the steps require a meta-data file associated with your character matrix.
The pipeline can be run without the meta-data file and still generate an output matrix which can be merged with Open Tree, but it won't distinguish between
asserted and inferred taxa states and will not print the literature sources for conflicting states. Your Phenoscape matrix should contain only one character
at a time.
'''

# First step of the pipeline is to convert the nexml matrix downloaded from Phenoscape to a tab-delimited version

sys.path.insert(0, os.path.join("tree_integration_pipeline", "1.nexml_to_txt_converter"))
import nexml_to_txt_converter

# Occasionally, there can be taxa with missing data/states (?). The second step is to remove them from the tab-delimited matrix
# After this step matrix pre-processing is complete

sys.path.insert(0, os.path.join("tree_integration_pipeline", "2.matrixpreprocessor"))
import matrixpreprocessor

# Third step is to remove conflict states '0&1' from the higher-level taxons (except species) from the pre-processed matrix
# These conflicts are due to conflicting statements by authors on the same character for same higher-level taxa

sys.path.insert(0, os.path.join("tree_integration_pipeline", "3.conflict_state_remover"))
import conflict_state_remover

# requesting users whether they have metadata or not
print '\n'
print 'requesting meta-data..'



while True:
    metainfo = raw_input('do you have meta deta? (Type y if you do, n if you do not):').lower()
# If the user has a separate meta-data file, the following streps within the if statement will be run
    if metainfo == 'y':
        print ('please move your metadata file to the folder named inputs')
        while True:
            metafile = raw_input('please type the full name (with extension) of the metadata xml file:')

            # Using error handling to prompt the user to correctly enter the file name
            try:
                p = open(os.path.join("inputs", metafile), 'r')
                p.close()
                break

            except (IOError):
                print "The file name you entered is wrong. It is not found in the inputs folder. Please retry!"
                continue

        # There is a part of the 3rd step which needs the meta-data file
        # This section will generate a file in statistics folder which provides literature sources for the conflicts found during the 3rd step
        # Therefore, user can find the reason for each conflicting state
        conflict_state_remover.conflict_state_remove(metafile)

        #The 4th step is to distinguish between inferred vs asserted states and replace the inferred presence with 2 instead of 1 and inferred absence with 3 instead of 0
        #This step will be run only if the meta-data is available

        sys.path.insert(0, os.path.join("tree_integration_pipeline", "4.inferred_state_replace"))
        import inferred_state_replace

        inferred_state_replace.inferred_state_replacer(metafile)

        break # when you get the correct input, break the while loop

    elif metainfo == 'n':
        print 'You do not have a meta-data file!!'
        print 'Without meta-data it is impossible to distinguish inferred vs asserted states, hence, step 4 will be skipped'
        print 'However the remaining steps will be continued'

        break # when you get the correct input, break the while loop

    # If the user provides a wrong input, prompt the user until getting the correct input
    else:
        print '\n'
        print 'incorrect input!!!!!!  please type y or n'
        #quit()
        continue # return back to the top of the while loop



# The 5th step is propagation; This step will propagate the character states of higher-level taxa to their species as higher level taxa cannot be
# mapped to open tree phylogeny

sys.path.insert(0, os.path.join("tree_integration_pipeline", "5.propagation"))
import propagation
propagation.propagation(metainfo)

# The final step of the pipeline is taxonomic reconciliation
# The Phenoscape supermatrix contains taxa based on Vertebrate Taxonomy  Ontology and the Open Tree taxa are based on NCBI
# The final steps matches the taxa between two sources and trasfers the data from VTO names to Open Tree names so the final matrix can be merged with Open Tree
# For this step, you need to input the name of the tree file downloaded from Open Tree

sys.path.insert(0, os.path.join("tree_integration_pipeline", "6.taxonomic_reconciliation"))
import taxonomic_reconciliation
taxonomic_reconciliation.taxonomic_reconciliation(metainfo)


print '\n'
print 'The pipeline is successfully completed. The final output matrices can be found in the outputs folder'
print 'These matrices contain Open Tree names and can be merged with the Open Tree phylogeny using a software like Mesquite'
print 'Moreover, the pipeline will generate a pre-processed version of your input Open Tree phylogeny (removing OT ids from the name)'
print 'This file is also in the outputs folder. Use this file for merging. Otherwise, merging will not work!!!'

