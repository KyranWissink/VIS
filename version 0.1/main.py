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
import functions as fc



####################################
#             Variables            #
####################################

# First argument is the input file and output file is name.csv
infile = sys.argv[1]
outfile = "output3.csv"



####################################
#          Read input file         #
####################################

with open(infile) as ifs:
    content = ifs.readlines()
content = [variant.strip() for variant in content]

    
    
####################################
#         Predict and write        #
####################################
    
def run(content):
    # Loop through all variants
    try:
        with open(outfile, "w") as ofs:
            ofs.write("Variant, Gene, 80-100, 50-80, 20-50, 5-20\n") # .csv headers
            for variant in content:
                try:
                    # Get the gene and transcript predictions
                    gene, predictions = fc.getResults(variant)
                    # Output it in csv with comma as delimiter
                    output = variant + "," + gene + "," + predictions["80"] + "," \
                        + predictions["50"] + "," + predictions["20"] + "," + \
                            predictions["5"] + "\n"
                    ofs.write(output)
                except:
                    # If anything fails, print the variant
                    print(variant + " failed!")
                    pass
                print("\n")
    except:
        print("Unable to open output file: " + outfile)
        sys.exit()


if __name__ == '__main__':
    run(content)





