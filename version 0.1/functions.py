# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:39:19 2021

@author: Kyran Wissink

Auxillary functions for the HGVSpredict.py
"""

####################################
#              Modules             #
####################################

import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import re
import requests, sys
from pyfaidx import Fasta
from hgvs.exceptions import HGVSInvalidVariantError


# Initialise HGVS to hg38 server
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()
am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name='GRCh38', 
                                        alt_aln_method='splign', 
                                        replace_reference=True)

# Reads the hg38 as a fasta
genome = Fasta('hg38.fa')

# Hold a small database of used gene to save time on ENSEMBL api
gene_list = []
gene_db = []


####################################
#             Functions            #
####################################

# The main control function that passes variants through other functions
def getResults(variant):
    print(variant)
    # Parse from HGVS to hg38
    try:
        genomic_variant, pos = HGVStohg38(variant)
    except:
        raise HGVSInvalidVariantError("not a valid HGVS variant.")
        return
    print(genomic_variant)
    
    # Get scores and gene from spliceAI
    try:
        gene, score, mrnapos = spliceAI(genomic_variant)
    except RuntimeError:
        print("SpliceAI was unable to handle " + variant + ".")
        return 
    print(score)
    
    # Find the location of the mutation in the transcripts
    try:
        location = locationFinder(genomic_variant, pos, gene)
    except RuntimeError:
        print("Location within gene not found. ")
        return
    print(location)
    
    # Make a dict of the delta score predictions
    predictions = {
        "80": "", 
        "50": "", 
        "20": "",
        "5": ""
    }
    
    # If any delta score is greater than 0.8
    if sum(i >= 0.8 for i in score.values()):
        try:
            predictions["80"] = predicter(score, mrnapos, location, 0.8)
        except RuntimeError:
            print("There was an error with the prediction of variant " + variant + ".")

    
    # If any delta score is between 0.5 and 0.8
    if sum(i >= 0.5 and i < 0.8 for i in score.values()):
        try:
            predictions["50"] = predicter(score, mrnapos, location, 0.5)
        except RuntimeError:
            print("There was an error with the prediction of variant " + variant + ".")

    
    # If any delta score is between 0.2 and 0.5
    if sum(i >= 0.2 and i < 0.5 for i in score.values()):
        try:
            predictions["20"] = predicter(score, mrnapos, location, 0.2)
        except RuntimeError:
            print("There was an error with the prediction of variant " + variant + ".")

        
    # If any delta score is between 0.05 and 0.2
    if sum(i >= 0.05 and i < 0.2 for i in score.values()):
        try:
            predictions["5"] = predicter(score, mrnapos, location, 0.05)
        except RuntimeError:
            print("There was an error with the prediction of variant " + variant + ".")

    print(predictions)
    # Return the data
    return gene, predictions
    


# Converts the HGVS cDNA variant to genomic hg38   
def HGVStohg38(hgvs_variant):
    var_c = hp.parse_hgvs_variant(hgvs_variant)
    
    # Use a previous version if the supplied version is not found
    # Loops through counting down until either found or not
    original_version = int(var_c.ac[-1])
    while True:
        try:
            var_g = am.c_to_g(var_c) # tries to convert cDNA to genomic
        except:
            # Try previous version if this one does not work
            version = int(var_c.ac[-1]) -1
            
            # If 9 other versions were tested, skip
            if version == original_version:
                print ("Variant " + hgvs_variant + " could not be parsed. Variant skipped.\n")
                break
            
            # Try later versions if version 0 is reached
            if version == 0:
                version = 9

            var_c.ac = var_c.ac[0:-1] + str(version) # Updates version
            continue
        break
    
    # Get the chromosome from NC number
    if str(var_g.ac[7]) != '0': # double digit chromosome
        chrom = str(var_g.ac[7]) + str(var_g.ac[8]) 
    else: # single digit chromosome
        chrom = str(var_g.ac[8])
    
    # X is chromosome 23
    if chrom == '23':
        chrom = 'X'
    
    # dels
    if "del" in hgvs_variant:
        # UNFINISHED #
        start = var_g.posedit.pos.start.base - 2 # -1 for extra nt, -1 for array start
        end = var_g.posedit.pos.end.base
        pos = str(start + 1)
        ref = str(genome["chr" + chrom][start:end]) # deleted sequence
        alt = ref[0] # first nucleotide of ref
    
    # dups
    elif "dup" in hgvs_variant:
        return
    
    #inserts
    elif "ins" in hgvs_variant:
        return
    
    # SNVs
    else:
        pos = str(var_g.posedit.pos.start.base)
        ref = str(var_g.posedit.edit.ref)
        alt = str(var_g.posedit.edit.alt)

    # Final genomic variant construction
    genomic_variant = chrom + "-" + pos + "-" + ref + "-" + alt
    return genomic_variant, pos # pos is later used to locate variant
    
    
    
# Runs the genomic variant through the spliceAI API
def spliceAI(variant):
    
    # Define the url
    url = 'https://spliceailookup-api.broadinstitute.org/spliceai/?hg=38&variant=' + variant
    
    # Pull data from API and format it into a list
    with requests.get(url) as r:
        data = r.json()["scores"][0].split("|")
    
    # Define the gene -- trim if necessary
    gene = data[0]
    if "ENSG" in gene:
        gene = re.search("([\w\d]+)\-", gene).group(1)
    
    # Make a dict of the scores
    score = {
        "acceptor_gain": float(data[1]),
        "acceptor_loss": float(data[2]),
        "donor_gain": float(data[3]),
        "donor_loss": float(data[4])
        }
    
    # Make a dict of the mRNA positions of all the scores
    mrnapos = {
        "acceptor_gain": int(data[5]),
        "acceptor_loss": int(data[6]),
        "donor_gain": int(data[7]),
        "donor_loss": int(data[8])
        }
    
    return gene, score, mrnapos



# Finds the location of the mutation (e.g. inside intron 11)
# Also adds the distance to the start and end of that intron/exon from the mut     
def locationFinder(variant, pos, gene):
    
    # If the gene was already pulled from the ENSEMBL database
    if gene in gene_list:
        index = gene_list.index(gene)
        decoded = gene_db[index]
    
    else:
        # Define server and query
        server = "http://rest.ensembl.org"
        ext = "/lookup/symbol/homo_sapiens/" + gene +  "?expand=1"
        r = requests.get(server+ext, headers={"Content-Type" : "application/json"})
        
        # Outputs an error and quits if the API is down
        if not r.ok:
          r.raise_for_status()
          sys.exit()
        
        # Output of ENSEMBL is a big json file
        decoded = r.json()
        
        # Add gene to local database
        gene_list.append(gene)
        gene_db.append(decoded)

        
    exoncount = 0 # counts the exons
    pos = int(pos)
    
    # Determine which transcript to use based on exon count
    transcript_list = []
    for i in range(0,len(decoded["Transcript"])):
        transcript_list.append(len(decoded["Transcript"][i]["Exon"]))
    TN = transcript_list.index(max(transcript_list))
    
    # Loop through the exons and find which boundary it is
    # Some genes are stored '3 to '5 and others '5 to '3 
    while True:
        
        exon_start = decoded["Transcript"][TN]["Exon"][exoncount]["start"]
        exon_end = decoded["Transcript"][TN]["Exon"][exoncount]["end"]
        next_exon_start = decoded["Transcript"][TN]["Exon"][exoncount + 1]["start"]
        # Within exon 5' to 3'
        if exon_start <= pos <= exon_end:
            location = {
                "exon": exoncount + 1, 
                "start": pos - exon_start, 
                "end": exon_end - pos
                }
            break
        
        # Within exon 3' to 5'
        elif exon_end <= pos <= exon_start:
            location = {
                "exon": exoncount + 1, 
                "start": pos - exon_end, 
                "end": exon_start - pos
                }
            break
            
        # Within intron 5' to 3'
        elif exon_end <= pos <= next_exon_start:
            location = {
                "intron": exoncount + 1,
                "start": pos - exon_end,
                "end": next_exon_start - pos
                }
            break
        
        # Within intron 3' to 5'
        elif next_exon_start <= pos <= exon_end:
            location = {
                "intron": exoncount + 1,
                "start": pos - next_exon_start,
                "end": exon_end - pos
                }
            break
        
        # If not found within exon or intron
        else:
            exoncount = exoncount + 1 # update exon count
            continue


    # Returns which exon or intron the mutation is inside of
    return location



# Predicts the effect of the mutation based on scores and location
def predicter(score, mrnapos, location, threshold):
    
    # Make it a bit simpler to work with
    acc_gain = score["acceptor_gain"]
    acc_loss = score["acceptor_loss"]
    don_gain = score["donor_gain"]
    don_loss = score["donor_loss"]

    
    ####################################
    #           Acceptor gain          #
    ####################################
    if acc_gain >= threshold:
        
        # Acceptor gain alone
        if acc_loss <= threshold and don_gain <= threshold and don_loss <= threshold:
            
            # If the mutation is in an exon
            if "exon" in location:
                bp = mrnapos["acceptor_gain"] + location["start"] # acceptor sites are at the beginning of introns
                
                # If the bp of the acceptor gain is before the start of the exon (negative because '5 )
                # That would lead to a partial insert
                if bp < 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(abs(bp)) + " bp ins"
                
                # Acceptor gain inside the exon would result in a partial del
                if bp > 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(bp) + " bp del"
                    
            # If the mutation is in an intron
            else:
                bp = mrnapos["acceptor_gain"] - location["end"] # acceptor sites are at the end of introns
                
                # if the acceptor gain is in an intron; partial insertion
                if bp < 0:
                    prediction = "Exon " + str(location["intron"] + 1) + " " + str(abs(bp)) + " bp ins"
                    
                # If the acceptor gain is in the exon; partial deletion
                if bp > 0:
                    prediction = "Exon " + str(location["intron"] + 1) + " " + str(bp) + " bp del"
                    
        # Acceptor gain + loss
        elif acc_loss >= threshold and don_gain <= threshold and don_loss <= threshold:

            # bp del or ins is based on how far acceptor gain is over loss
            # negative is upstream
            bp = mrnapos["acceptor_gain"] - mrnapos["acceptor_loss"]

            # If the mutation is in an exon
            if "exon" in location:

                # upstream acceptor gain over loss means an insert
                if bp < 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(abs(bp)) + " bp ins"
                    
                # Downstream acceptor gain over loss means del
                if bp > 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(bp) + " bp del"
                    
            # if the mutation is in an intron
            else:

                # upstream acceptor gain in intron extends the downstream exon
                if bp < 0:
                    prediction = "Exon " + str(location["intron"] + 1) + " " + str(abs(bp)) + " bp ins"
                # Downstream acceptor gain in intron dels the downstream exon
                if bp > 0:
                    prediction = "Exon " + str(location["intron"] + 1) + " " + str(bp) + " bp del"
                    
        # Acceptor gain + loss + donor loss/gain
        else:
            prediction = "complex"
                    
                    
                    
    ####################################
    #           Acceptor loss          #
    ####################################
    elif acc_loss >= threshold:
        
        # Acceptor loss alone
        if don_gain <= threshold and don_loss <= threshold:
            
            # If the mutation is in an exon that would be skipped
            if "exon" in location:
                prediction = "Exon " + str(location["exon"]) + " skip"
                    
            # If the mutation is in an intron the exon after is skipped
            else:
                prediction = "Exon " + str(location["intron"] + 1) + " skip"
        
        # Acceptor loss + donor loss
        elif don_loss >= threshold and don_gain <= threshold:
            
            # Exon skip if donor loss is downstream (positive 3') of acceptor loss
            if mrnapos["donor_loss"] > mrnapos["acceptor_loss"]:
                
                if "exon" in location:
                    prediction = "Exon " + str(location["exon"]) + " skip"
                    
                else:
                    prediction = "Exon " + str(location["intron"]) + " skip"
            
            # Double exon skip or intron retention if donor loss is upstream (negative 5') of acceptor loss
            elif mrnapos["donor_loss"] < mrnapos["acceptor_loss"]:
                
                if "exon" in location:
                    prediction = "Exon " + str(location["exon"] - 1) + " & " + str(location["exon"])  + " skip" + \
                                " OR intron " + str(location["exon"] - 1 + " retention")
                
                else:
                    prediction = "Exon " + str(location["intron"] - 1) + " & " + str(location["intron"])  + " skip" + \
                                " OR intron " + str(location["intron"] - 1 + " retention")
                                
            # Exon skip if donor loss is downstream (positive 3')
            else:
                
                if "exon" in location:
                    prediction = "Exon " + str(location["exon"]) + " skip"
                else:
                    prediction = "Exon " + str(location["intron"]) + " skip"
        
        # Acceptor loss + donor loss + donor gain
        else:
            prediction = "complex"


    ####################################
    #             Donor gain           #
    ####################################
    elif don_gain >= threshold:
        
        # Donor gain alone
        if don_loss <= threshold:
            
            # If the mutation is in an exon
            if "exon" in location:
                bp = mrnapos["donor_gain"] - location["end"] # donor sites are at the end of introns
                
                # If the bp of the donor gain is before the end of the exon (negative because '5 )
                # That would lead to a partial delete
                if bp < 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(abs(bp)) + " bp del"
                
                # Donor gain inside the intron would result in a partial del
                if bp > 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(bp) + " bp ins"
                    
            # If the mutation is in an intron
            else:
                bp = abs(mrnapos["donor_gain"]) - abs(location["start"]) # donor sites are at the beginning of introns
                
                # if the donor gain is in an intron
                if bp < 0:
                    prediction = "Exon " + str(location["intron"]) + " " + str(abs(bp)) + " bp ins"
                    
                # If the acceptor gain is in the exon
                if bp > 0:
                    prediction = "Exon " + str(location["intron"]) + " " + str(bp) + " bp del"
            
        # Donor gain + loss            
        elif don_loss >= threshold:
            
            # bp del or ins is based on how far donor gain is over loss
            # negative is upstream
            bp = mrnapos["donor_gain"] - mrnapos["donor_loss"]
            
            # If the mutation is in an exon
            if "exon" in location:
                
                # upstream donor gain over loss means del
                if bp < 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(abs(bp)) + " bp del"
                    
                # Downstream donor gain over loss means an insert
                elif bp > 0:
                    prediction = "Exon " + str(location["exon"]) + " " + str(bp) + " bp ins"
                    
            # if the mutation is in an intron
            else:
                
                # upstream donor gain in intron dels the downstream exon
                if bp < 0:
                    prediction = "Exon " + str(location["intron"]) + " " + str(abs(bp)) + " bp del"
                    
                # Downstream acceptor gain in intron extends the downstream exon
                elif bp > 0:
                    prediction = "Exon " + str(location["intron"]) + " " + str(bp) + " bp ins"
        
            
            
    
    ####################################
    #             Donor loss           #
    ####################################
    elif don_loss >= threshold:
        # Donor loss alone
        # If the mutation is in an exon that would be skipped
        if "exon" in location:
            prediction = "Exon " + str(location["exon"]) + " skip"
                
        # If the mutation is in an intron the exon before is skipped
        else:
            prediction = "Exon " + str(location["intron"]) + " skip"
    
    return prediction
                       
        
        

    







