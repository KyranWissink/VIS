# -*- coding: utf-8 -*-
"""
Created on Sat Mar 13 12:18:03 2021

@author: Kyran Wissink

Input: any file where HGVS mutations are seperated by newlines.

Uses the spliceAI API.

Installed dependencies:
    hgvs
    pyfaidx
    pyensemblrest
    pandas

Output: A .csv file with predictions of the effects of the mutations
        on the transcript.
"""

####################################
#              Modules             #
####################################

import sys, os
import HGVSfunctions as hf
import pandas as pd
import argparse



####################################
#          Argument setup          #
####################################

parser = argparse.ArgumentParser()
parser.add_argument('-i', "--input", help='input file of HGVS variants in plain text', required=True)
parser.add_argument('-o', "--output", help='output csv file of the variants', required=True)
parser.add_argument('-g', "--genome", help='genome version (hg19 or hg38)', required=True)
parser.add_argument('-p', "--pref_tr", help='OPTIONAL preffered transcripts txt')
args = parser.parse_args()



####################################
#         Predict and write        #
####################################
    
if __name__ == '__main__':
    
    print("\nStarting " + sys.argv[0] + ".")
    print("Initialising programme...")
    
    infile = args.input
    outfile = args.output
    
    with open(infile) as ifs:
        content = ifs.readlines()
    content = [variant.strip() for variant in content]
    
    if args.pref_tr:
        hf.validate_variants(content, args.pref_tr)
        
    hf.Predicter.grch = args.genome
    content = hf.Predicter.initialise(content)
    
    print("\rInitialisation complete.\n")
    print("Running variants...")
    
    write = [] # final output 
    for variant in content:
        print("Working on: " + variant + ".")
        output = hf.Predicter.run(variant)
        output = hf.clean(output) 
        write.append(output)
    
    write = pd.DataFrame(data=write)
    
    try:
        write.to_csv(outfile, index=False)
    except PermissionError:
        raise PermissionError("Unable to open the output file.\n" + 
                              "File is likely open in another programme.")
    
    print("\nSuccessfully processed " + str(len(write)) + " out of " + str(len(content)) + " variants.")
    print("Results can be found at " + os.path.abspath(os.getcwd()) + "/" + args.output)

