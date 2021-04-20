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


def check_arg(arg):
    """
    Parameters
    ----------
    arg : The argument's path file.

    Raises
    ------
    FileNotFoundError
        When the file is not found and exits the programme.

    """
    if not os.path.exists(arg):
        raise FileNotFoundError ("File not found for: %s." % arg)
        sys.exit()



####################################
#         Predict and write        #
####################################
    
if __name__ == '__main__':
    
    print("\nStarting " + sys.argv[0] + ".")
    print("Initialising programme...")


    # Check arguments
    check_arg(args.input)
    check_arg(args.output)
    check_arg(args.genome)
    if args.pref_tr:
        check_arg(args.pref_tr)
    
    
    # Read content and strip newline characters (\r\n)
    with open(args.input) as ifs:
        content = ifs.readlines()
    content = [variant.strip() for variant in content]
    
    
    # Validate if preferred transcripts are supplied
    if args.pref_tr:
        hf.validate_variants(content, args.pref_tr)
    
    
    # Initialise genome and content
    hf.Predicter.grch = args.genome
    content = hf.Predicter.initialise(content)
    
    print("Initialisation complete.\n")
    print("Running variants...")
    
    write = [] # final output 
    
    
    # Run every variant
    try:
        for variant in content:
            print("Working on: " + variant + ".")
            output = hf.Predicter.run(variant) # Predict
            output = hf.clean(output) # Do not include everything from the class
            write.append(output) # Add to final output
    except Exception:
        print("variant skipped: %s." % variant)
    
    
    # Convert results to pandas dataframe
    write = pd.DataFrame(data=write)
    
    
    # Write to .csv file with pandas
    try:
        write.to_csv(args.output, index=False)
    except PermissionError:
        raise PermissionError("Unable to open the output file.\n" + 
                              "File is likely open in another programme.")
        sys.exit()
    
    
    # Finalise and exit
    print("\nSuccessfully processed " + str(len(write)) + " out of " + str(len(content)) + " variants.")
    print("Results can be found at " + os.path.abspath(os.getcwd()) + "/" + args.output)

