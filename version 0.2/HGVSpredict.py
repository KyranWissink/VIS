# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 12:18:03 2021

@author: Kyran Wissink

Input: any file where HGVS mutations are seperated by newlines.

Uses spliceAI and ENSEMBL APIs. When either of these are down,
the script will not run

Output: A .csv file with predictions of the effects of the mutations
        on the transcript.
"""

####################################
#              Modules             #
####################################

import sys
import HGVSfunctions as hf



####################################
#            Initialise            #
####################################

# First argument is the input file and output file is name.csv
infile = sys.argv[1]
outfile = "output\output3.csv"



####################################
#          Read input file         #
####################################

with open(infile) as ifs:
    content = ifs.readlines()
content = [variant.strip() for variant in content]

genome = sys.argv[2]
    
content, genome, grch = hf.initialise(content, genome)



####################################
#         Predict and write        #
####################################
    
def run(content):
    with open(outfile, "w") as ofs:
        # .csv headers
        ofs.write("Variant,Gene,80-100,50-80,20-50,5-20,Warning\n") 
        
        # Loop through all variants
        for variant in content:
                print(variant)
                # Get the gene and transcript predictions
                output = hf.predicter(variant, genome, grch)
                
                # Output it in csv with comma as delimiter
                results = variant + "," + output.gene + ","
                for prediction in output.predictions:
                    results += output.predictions[prediction] + ","
                    
                if hasattr(output, 'warn'):
                    results += output.warn
                results += "\n"
                
                ofs.write(results)


if __name__ == '__main__':
    run(content)





