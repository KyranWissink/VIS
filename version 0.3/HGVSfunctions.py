# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 18:03:18 2021

@author: Kyran Wissink

Auxillary functions for the HGVSpredict.py
"""
##########################################################################
#                           MODULES AND VARIABLES                        #
##########################################################################

import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import requests, re
import sys, time
from pyfaidx import Fasta
from hgvs.exceptions import HGVSInvalidVariantError
from functools import lru_cache
import hgvs.exceptions as he

# HGVS server 
hdp = hgvs.dataproviders.uta.connect()
hp = hgvs.parser.Parser()



##########################################################################
#                               INITIALISATION                           #
##########################################################################

def initialise(content, genome):
    genome, grch = getGenome(genome)
    content = validateVariants(content)
    
    return content, genome, grch

# Get preferred transcripts for all variants
def validateVariants(content):
    
    try:
        with open ("input/preferred_transcript.txt") as ifs:
            preferred = ifs.readlines()
        preferred = [variant.strip() for variant in preferred]
        
        # Verify transcripts
        for index in range(0,len(content)):
            match = re.search("(\w+)\.\d+(.*)", content[index])
            nm = match.group(1)
            rest = match.group(2)
            transcript = [x for x in preferred if nm in x]
            if len(transcript) == 1:
                transcript = transcript[0]
                content[index] = transcript + rest
            else:
                transcript = content[index]
                
    # Just return content when there are no preferred transcripts
    except FileNotFoundError:
        return content
    return content


# Get the correct genome fasta
def getGenome(genome):
    if '19' in genome or '37' in genome:
        grch = "GRCh37"
        genome = Fasta('input/hg19.fa')
    elif '38' in genome:
        grch = "GRCh38"
        genome = Fasta('input/hg38.fa')
    else:
        raise RuntimeError("Invalid genome: " + genome + ".")
        sys.exit()
        
    return genome, grch
    


# Get the correct hgvs assembly mapper for genome version
def getAssemblyMapper(grch):
    am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=grch, 
                                            alt_aln_method='splign', 
                                            replace_reference=True)
    return am



##########################################################################
#                             GENOMIC VARIANT                            #
##########################################################################

# Convert HGVS to genomic notation
def getGenomic(var):
    am = getAssemblyMapper(var.grch)
    start = var.hgvs
    
    while True:    
        try:
            var_c = hp.parse_hgvs_variant(var.hgvs)
            var_g = am.c_to_g(var_c) # tries to convert cDNA to genomic
            var.chrom = re.search("([0-9]+)."  ,var_g.ac).group(1)[4:6] # Chrom from NC number
        
        # Try an older version of the NM if this version does not work
        except he.HGVSDataNotAvailableError:
            match = re.search("(^\w+\.\d+)(.*)", var.hgvs)
            NM = str(match.group(1))
            newversion = int(NM[-1]) + 1
            NM = NM[0:-1] + str(newversion)[-1]
            var.hgvs = NM + str(match.group(2))
            
            # Make it fail if it looped ten times
            if var.hgvs == start: 
                return
            var.warn = "HGVS error: genomic variant conversion was manually performed"
            continue
        
        # If all fails
        except Exception:
            return
        
        break

    # Remove zeros from one digit chromosomes
    if var.chrom[0] == '0':
        var.chrom = var.chrom[1]
        
    if var.chrom == '23':
        var.chrom = 'X'
    
    # dels
    if "del" in var.hgvs:
        delGenomic(var, var_g)
        
    # dups
    elif "dup" in var.hgvs:
        dupGenomic(var, var_g)
        
    # ins
    elif "ins" in var.hgvs:
        insGenomic(var, var_g)
        
    #SNV
    else:
        snvGenomic(var, var_g)
        
# Use the genome fasta to determine ref and alt for dels
def delGenomic(var, var_g):
    start = var_g.posedit.pos.start.base - 2
    end = var_g.posedit.pos.end.base - 1
    ref = str(var.genome["chr" + var.chrom]\
        [ start - 1: end ])
    alt = ref[0]
    
    var.pos = str(start)
    var.genomic = var.chrom + "-" + var.pos + "-" + ref + "-" + alt
    
# Use the genome fasta to determine ref and alt for dups
# UNFINISHED
def dupGenomic(var, var_g):
    start = var_g.posedit.pos.start.base
    end = var_g.posedit.pos.end.base
    
        
    dup = str(var.genome["chr" + var.chrom]\
                  [var.pos - 1: end])
    ref = dup[0]
    alt = dup[0]+dup
    
    var.pos = str(start)
    var.genomic = var.chrom + "-" + var.pos + "-" + ref + "-" + alt
        
    
# UNFINISHED
def insGenomic(var):
    var.genomic = ""

# Get genomic variant from the hgvs module results
def snvGenomic(var, var_g):
    var.pos = str(var_g.posedit.pos.start.base)
    var.ref = str(var_g.posedit.edit.ref)
    var.alt = str(var_g.posedit.edit.alt)
    
    var.genomic = var.chrom + "-" + var.pos + "-" + var.ref + "-" + var.alt



##########################################################################
#                               TRANSCRIPTS                              #
##########################################################################  
        

# Gets all transcripts of the gene from ENSEMBL using the REST API
def getEnsembl(var):
    # Define server and query
    if "38" in var.grch:
        server = "http://rest.ensembl.org"
        ext = "/lookup/symbol/homo_sapiens/" + var.gene +  "?expand=1"
        transcripts = requests.get(server+ext, headers={"Content-Type" : "application/json"}, timeout=10)
    else:
        server = "http://grch37.rest.ensembl.org"
        ext = "/lookup/symbol/homo_sapiens/" + var.gene + "?content-type=application/json;expand=1"
        transcripts = requests.get(server+ext, timeout=10)
        

    
    # Outputs an error and quits if the API is down
    if not transcripts.ok:
      transcripts.raise_for_status()
      sys.exit()
    
    # Output of ENSEMBL is a big json file
    var.transcripts = transcripts.json()["Transcript"]
    


##########################################################################
#                               LOCATIONS                                #
##########################################################################     

# Sketchy code to get all locations of the mutation in valid transcripts
def getLocation(var):
    location_counter = {}
    location_list = []
    for transcript in var.transcripts:
        location = locationFinder(transcript, var.pos)
        if location == "ERROR":
            continue
        if location['number'] not in location_counter:
            location_list.append(location)
            location_counter[location['number']] = 1
        else:
            location_counter[location['number']] += 1
            
    # Use most common location
    try:
        location_number = max(location_counter, key = location_counter.get)
    except ValueError:
        print("Location not found for " + var.hgvs)
    var.location = next(item for item in location_list if item["number"] == location_number)

    

def locationFinder(transcript, pos):
    
    # Loop through the exons and find which boundary it is
    # Some genes are stored '3 to '5 and others '5 to '3 
    exoncount = 0
    pos = int(pos)
    try:
        while True:
            exon_start = transcript["Exon"][exoncount]["start"]
            exon_end = transcript["Exon"][exoncount]["end"]
            next_exon_start = transcript["Exon"][exoncount + 1]["start"]
            
            # Within exon 5' to 3'
            if exon_start <= pos <= exon_end:
                location = {
                    "region": "exon",
                    "number": exoncount + 1, 
                    "start": pos - exon_start, 
                    "end": exon_end - pos
                    }
                break
            
            # Within exon 3' to 5'
            elif exon_end <= pos <= exon_start:
                location = {
                    "region": "exon",
                    "number": exoncount + 1, 
                    "start": pos - exon_end, 
                    "end": exon_start - pos
                    }
                break
                
            # Within intron 5' to 3'
            elif exon_end <= pos <= next_exon_start:
                location = {
                    "region": "intron",
                    "number": exoncount + 1,
                    "start": pos - exon_end,
                    "end": next_exon_start - pos
                    }
                break
            
            # Within intron 3' to 5'
            elif next_exon_start <= pos <= exon_end:
                location = {
                    "region": "intron",
                    "number": exoncount + 1,
                    "start": pos - next_exon_start,
                    "end": exon_end - pos
                    }
                break
            
            # If not found within exon or intron
            else:
                exoncount = exoncount + 1 # update exon count
                continue
    except IndexError:
        location = "ERROR"

    # Returns which exon or intron the mutation is inside of
    return location
    

    
##########################################################################
#                                SPLICEAI                                #
##########################################################################      

def getSpliceAI(var):

    # Define the url
    if "37" in var.grch:
        spliceaigrch = "hg=37"
    else:
        spliceaigrch = "hg=38"
        
    url = 'https://spliceailookup-api.broadinstitute.org/spliceai/?' + spliceaigrch + '&variant=' + var.genomic
    
    # Pull data from API and format it into a list
    r = requests.get(url, timeout=30)
        
    if not r.ok:
        print("SpliceAI did not respond. Trying again in 10 seconds.")
        time.sleep(10)
        r = requests.get(url, timeout=30)
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    # Get score data and format to list
    data = r.json()["scores"][0].split("|")
    var.gene = data[0]
    
    # Trim genes, necessary for dels and dups 
    if "-" in var.gene:
        var.gene = re.search("^(\w+)--", var.gene).group(1)

    # Make a dict of the scores
    var.scores = {
        "acceptor_gain": float(data[1]),
        "acceptor_loss": float(data[2]),
        "donor_gain": float(data[3]),
        "donor_loss": float(data[4])
        }
    
    # Make a dict of the mRNA positions of all the scores
    var.mrnapos = {
        "acceptor_gain": int(data[5]),
        "acceptor_loss": int(data[6]),
        "donor_gain": int(data[7]),
        "donor_loss": int(data[8])
        }           



##########################################################################
#                              PREDICTIONS                               #
##########################################################################    

# Make predictions for all the delta score cutoffs
def getPredictions(var):
    var.predictions = {
        0.8: "",
        0.5: "",
        0.2: "",
        0.05: ""
        }
    
    # Loop through the cutoffs
    maxcut = 1.0
    for cutoff in var.predictions:
        if sum(i >= cutoff and i < maxcut for i in var.scores.values()):
            var.predictions[cutoff] = predict(var, cutoff)
        maxcut = cutoff
        
        
# Predict splicing effect based on spliceai scores and position
def predict(var, threshold):
    
    # Make it a bit simpler to work with
    acc_gain = var.scores["acceptor_gain"]
    acc_loss = var.scores["acceptor_loss"]
    don_gain = var.scores["donor_gain"]
    don_loss = var.scores["donor_loss"]
    
    
    
    ####################################
    #           Acceptor gain          #
    ####################################
    if acc_gain >= threshold:
        
        # Acceptor gain alone
        if acc_loss <= threshold and don_gain <= threshold and don_loss <= threshold:
            
            # If the mutation is in an exon
            if var.location["number"] == "exon":
                
                # acceptor sites are at the beginning of introns
                bp = var.mrnapos["acceptor_gain"] + var.location["start"] 
                
                # If the bp of the acceptor gain is before the start of the exon (negative because '5 )
                # That would lead to a partial insert
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp ins"
                
                # Acceptor gain inside the exon would result in a partial del
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp del"
                    
            # If the mutation is in an intron
            else:
                bp = var.mrnapos["acceptor_gain"] - var.location["end"] # acceptor sites are at the end of introns
                
                # if the acceptor gain is in an intron; partial insertion
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"] + 1) + " " + str(abs(bp)) + " bp ins"
                    
                # If the acceptor gain is in the exon; partial deletion
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"] + 1) + " " + str(bp) + " bp del"
                    
        # Acceptor gain + loss
        elif acc_loss >= threshold and don_gain <= threshold and don_loss <= threshold:

            # bp del or ins is based on how far acceptor gain is over loss
            # negative is upstream
            bp = var.mrnapos["acceptor_gain"] - var.mrnapos["acceptor_loss"]

            # If the mutation is in an exon
            if var.location["number"] == "exon":

                # upstream acceptor gain over loss means an insert
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp ins"
                    
                # Downstream acceptor gain over loss means del
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp del"
                    
            # if the mutation is in an intron
            else:

                # upstream acceptor gain in intron extends the downstream exon
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"] + 1) + " " + str(abs(bp)) + " bp ins"
                # Downstream acceptor gain in intron dels the downstream exon
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"] + 1) + " " + str(bp) + " bp del"
                    
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
            if var.location["number"] == "exon":
                prediction = "Exon " + str(var.location["number"]) + " skip"
                    
            # If the mutation is in an intron the exon after is skipped
            else:
                prediction = "Exon " + str(var.location["number"] + 1) + " skip"
        
        # Acceptor loss + donor loss
        elif don_loss >= threshold and don_gain <= threshold:
            
            # Exon skip if donor loss is downstream (positive 3') of acceptor loss
            if var.mrnapos["donor_loss"] > var.mrnapos["acceptor_loss"]:
                
                if var.location["number"] == "exon":
                    prediction = "Exon " + str(var.location["number"]) + " skip"
                    
                else:
                    prediction = "Exon " + str(var.location["number"]) + " skip"
            
            # Double exon skip or intron retention if donor loss is upstream (negative 5') of acceptor loss
            elif var.mrnapos["donor_loss"] < var.mrnapos["acceptor_loss"]:
                
                if var.location["number"] == "exon":
                    prediction = "Exon " + str(var.location["number"] - 1) + " & " + str(var.location["number"])  + " skip" + \
                                " OR intron " + str(var.location["number"] - 1 + " retention")
                
                else:
                    prediction = "Exon " + str(var.location["number"] - 1) + " & " + str(var.location["number"])  + " skip" + \
                                " OR intron " + str(var.location["number"] - 1 + " retention")
                                
            # Exon skip if donor loss is downstream (positive 3')
            else:
                
                if var.location["number"] == "exon":
                    prediction = "Exon " + str(var.location["number"]) + " skip"
                else:
                    prediction = "Exon " + str(var.location["number"]) + " skip"
        
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
            if var.location["number"] == "exon":
                
                # donor sites are at the end of introns
                bp = var.mrnapos["donor_gain"] - var.location["end"] 
                
                # If the bp of the donor gain is before the end of the exon (negative because '5 )
                # That would lead to a partial delete
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp del"
                
                # Donor gain inside the intron would result in a partial del
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp ins"
                    
            # If the mutation is in an intron
            else:
                
                # donor sites are at the beginning of introns
                bp = abs(var.mrnapos["donor_gain"]) - abs(var.location["start"]) 
                
                # if the donor gain is in an intron
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp ins"
                    
                # If the acceptor gain is in the exon
                if bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp del"
            
        # Donor gain + loss            
        elif don_loss >= threshold:
            
            # bp del or ins is based on how far donor gain is over loss
            # negative is upstream
            bp = var.mrnapos["donor_gain"] - var.mrnapos["donor_loss"]
            
            # If the mutation is in an exon
            if var.location["number"] == "exon":
                
                # upstream donor gain over loss means del
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp del"
                    
                # Downstream donor gain over loss means an insert
                elif bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp ins"
                    
            # if the mutation is in an intron
            else:
                
                # upstream donor gain in intron dels the downstream exon
                if bp < 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(abs(bp)) + " bp del"
                    
                # Downstream acceptor gain in intron extends the downstream exon
                elif bp > 0:
                    prediction = "Exon " + str(var.location["number"]) + " " + str(bp) + " bp ins"
        
        
            
    ####################################
    #             Donor loss           #
    ####################################
    elif don_loss >= threshold:
        # Donor loss alone
        # If the mutation is in an exon that would be skipped
        if var.location["number"] == "exon":
            prediction = "Exon " + str(var.location["number"]) + " skip"
                
        # If the mutation is in an intron the exon before is skipped
        else:
            prediction = "Exon " + str(var.location["number"]) + " skip"
    
    return prediction



##########################################################################
#                             CLASS DEFINITION                           #
##########################################################################      


class predicter(object):
    def __init__(self, hgvs, genome, grch):
        self.hgvs = hgvs
        self.genome = genome
        self.grch = grch
        getGenomic(self)
        getSpliceAI(self)
        getEnsembl(self)
        #validateTranscripts(self)
        getLocation(self)
        getPredictions(self)

    
    def __tuple__(self):
        return tuple(self)
        
