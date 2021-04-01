# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 12:18:03 2021

@author: Kyran Wissink

Input: any file where HGVS mutations are seperated by newlines.

Uses spliceAI and ENSEMBL APIs.

Output: A .csv file with predictions of the effects of the mutations
        on the transcript.
"""

####################################
#              Modules             #
####################################

import sys
import HGVSfunctions as hf


####################################
#         Predict and write        #
####################################
    
def run(content):
    
    with open(outfile, "w") as ofs:
        # .csv headers
        ofs.write("Variant,Gene,80-100,50-80,20-50,5-20,Warning\n") 
        
        # Loop through all variants
        for variant in content:
            
            try:
                print(variant)
                
                # Get the gene and transcript predictions
                output = hf.predicter(variant, genome, grch)
                
                print(output.genomic)
                print(output.predictions)
                print("\n")
                # Output it in csv with comma as delimiter
                results = variant + "," + output.gene + ","
                for prediction in output.predictions:
                    results += output.predictions[prediction] + ","
                
                # Write warnings if there are any
                if hasattr(output, 'warn'):
                    results += output.warn
                results += "\n"
                
                ofs.write(results)

            except Exception:
                ofs.write(variant + "," + "ERROR" + "\n")


if __name__ == '__main__':
        
    # First argument is the input file and output file is name.csv
    infile = sys.argv[1]
    outfile = "output/output3.csv"
    
    with open(infile) as ifs:
        content = ifs.readlines()
    content = [variant.strip() for variant in content]
    
    genome = sys.argv[2]
    content, genome, grch = hf.initialise(content, genome)
    
    run(content)





