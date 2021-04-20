# -*- coding: utf-8 -*-
"""
Created on Wed Apr 7 11:38:32 2021

@author: Kyran Wissink

Auxiliary functions to HGVSpredict.py
"""
##########################################################################
#                                  MODULES                               #
##########################################################################

import hgvs.dataproviders.uta
import hgvs.assemblymapper
import hgvs.parser
import requests
import re, sys, time
from pyfaidx import Fasta
import hgvs.exceptions as he
from ensemblrest import EnsemblRest
ensRest = EnsemblRest()



##########################################################################
#                                  METHODS                               #
##########################################################################

def validate_variants(content):
    """
    Parameters
    ----------
    content : List of all variants in hgvs format, parsed from the user 
            input .txt file in the command line options.

    Returns
    -------
    content : List of all variants in hgvs format, now validated with 
            the preferred transcripts.

    """

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
            
    return content



def clean(var):
    """
    Function
    --------
    Makes it so that unimportant information such as the VEP and ENSEMBL 
    transcripts do not go in the output file.
    
    Parameters
    ----------
    var : The variant most recently processed.

    Returns
    -------
    output : All relevant information from the variant to be parsed to 
        the final csv file by pandas. 

    """
    
    output = {}
    output["HGVS"] = var["hgvs"]
    output["Gene"] = var["gene"]
    output["Chrom"] = var["chrom"]
    output["Genomic"] = var["genomic"]
    output["80-100%"] = var["prediction"][0.8]
    output["50-80%"] = var["prediction"][0.5]
    output["20-50%"] = var["prediction"][0.2]
    output["5-20%"] = var["prediction"][0.05]
    
    return output



class Predicter():
    
    ##########################################################################
    #                               INITIALISATION                           #
    ########################################################################## 
    
    @staticmethod
    def get_genome(grch):
        """
        Parameters
        ----------
        grch : The genome version supplied by the user from the command line.

        Raises
        ------
        RuntimeError
            When the genome is not either 37 or 38.

        Returns
        -------
        genome : The full genome in fasta format, seperated by chromosomes.
        grch : Standardised genome version.

        """
        
        if '19' in grch or '37' in grch:
            grch = "GRCh37"
            genome = Fasta('input/hg19.fa')
            
        elif '38' in grch:
            grch = "GRCh38"
            genome = Fasta('input/hg38.fa')
            
        else:
            raise RuntimeError("Invalid genome: " + grch + ".")
            sys.exit()
            
        return genome, grch
    
    
    
    @staticmethod
    def validate_variants(content, pref_tr):
        """
        Function
        --------
        Validates the transcripts by using the preferred HGVS transcript 
        variant if supplied. 
        
        This method is only called when the user inputs a preferred transcript
        file in the command line options.
        
        Parameters
        ----------
        content : List of all variants in hgvs format.

        Returns
        -------
        content : List of all variants in hgvs format, now validated with 
                the preferred transcripts.

        """

        with open (pref_tr) as ifs:
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
                
        return content
    
    

    @classmethod
    def get_hgvs_assembly(cls):
        """
        Function
        --------
        To set up the HGVS server and correct assembly based on the human 
        genome version. Gets the HG version from the class attribute.

        """
        
        hdp = hgvs.dataproviders.uta.connect()
        cls.hp = hgvs.parser.Parser()
        cls.am = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name=cls.grch, 
                                            alt_aln_method='splign', 
                                            replace_reference=True)
        
        
        
    @classmethod
    def initialise(cls, content):
        """
        Function
        ----------
        Initialises the class prior to running any variants. It sets the 
        genome Fasta as a class variable, sets the GRCh in the right format 
        as a class variable, and validates the variants.
        
        """
        
        cls.genome, cls.grch = cls.get_genome(cls.grch)
        cls.get_hgvs_assembly()
        
        return content
        
        
        
    ##########################################################################
    #                               PER-VARIANT RUN                          #
    ########################################################################## 
    
    @classmethod
    def run(cls, variant):
        """
        Function
        -------
        Called from external script to run a variant through all methods.
        Gets the VEP and ENSEMBL transcript information;
        Converts the HGVS variant to a genomic variant;
        Gets the SPLICEAI scores and positions;
        Predicts the effect of the mutation based on the above.

        Returns
        -------
        The variant dict with all the above information.

        """
        
        var = {}
        var["hgvs"] = variant
        
        try:
            var["vep"] = ensRest.getVariantConsequencesByHGVSnotation\
                (species="human", hgvs_notation=var["hgvs"])[0]
        except Exception as e:
            print(e)
            return var
        
        try:
            var["gene"] = cls.get_gene(var["vep"]["transcript_consequences"])
            var["chrom"] = var["vep"]["seq_region_name"]
        except Exception as e:
            print(e)
            return var
         
        try:
            var["lookup"] = ensRest.getLookupBySymbol\
                (species="homo_sapiens", symbol=var["gene"], expand=1)
        except Exception as e:
            print(e)
            return var
        
        try:
            var = cls.get_genomic(var)
        except Exception as e:
            print(e)
            return var
        
        try:
            var["scores"], var["mrnapos"] = cls.get_spliceai(var["genomic"], cls.grch)
        except Exception as e:
            print(e)
            return var
            
        try:
            var["location"] = cls.get_location(var["lookup"]["Transcript"], var["pos"])
        except Exception as e:
            print(e)
            return var
        
        try:
            var["prediction"] = cls.get_predictions(var)
        except Exception as e:
            print(e)
            return var
    
        return var
            
            
    
    @staticmethod
    def get_gene(vep):
        """
        Function
        --------
        Gets the gene from the VEP transcripts. 
        Achieves this by comparing all the genes that VEP outputs for this 
        variant and then proceeds to take the most common one.
        
        Parameters
        ----------
        vep : The ENSEMBL VEP["Transcript"] of the variant

        Returns
        -------
        gene : The gene symbol of the variant.

        """
        
        transcripts = len(vep)
        genes = {}
        
        for x in range(0, transcripts):
            if 'distance' in vep[x]: # Exclude genes that are far away
                continue
            gene = vep[x]["gene_symbol"]
            if gene in genes:
                genes[gene] += 1
            else:
                genes[gene] = 1
                
        gene = max(genes, key=genes.get)
        return gene
    
    
    
    ##########################################################################
    #                              GENOMIC VARIANT                           #
    ##########################################################################    
    
    @classmethod
    def get_genomic(cls, var):
        """
        Function
        --------
        Serves as a controlling method to get the genomic variant from the 
        cDNA HGVS variant. Calls upon different methods based on dels, dups,
        ins, or SNVs. 

        Returns
        -------
        var : updated variant, with pos and genomic now.

        """
    
        # dels
        if "del" in var["hgvs"]:
            var["pos"], ref, alt = cls.del_genomic(var, cls.genome)
            
        # dups
        elif "dup" in var["hgvs"]:
            var["pos"], ref, alt = cls.dup_genomic(var, cls.genome)
            
        # ins
        elif "ins" in var["hgvs"]:
            var["pos"], ref, alt = cls.ins_genomic(var)
            
        #SNV
        else:
            var["pos"], ref, alt = cls.snv_genomic(var)
            
        var["genomic"] = var["chrom"] + "-" + var["pos"] + "-" + ref + "-" + alt
            
        return var
    
    
    
    @staticmethod
    def dup_genomic(var, genome):
        """
        Function
        --------
        Gets the genomic variant from the HGVS variant for dup variants.
        
        Parameters
        ----------
        var : Variant dict.
        genome : The genome annotation Fasta.

        Returns
        -------
        pos : Genomic position of the mutation.
        ref : Ref at that genomic position.
        alt : Alt that has occured due to the mutation.

        """
        positions = [var["vep"]["start"], var["vep"]["end"]]
        start = min(positions)
        end = max(positions)
        ref = str(genome["chr" + str(var["chrom"])]\
            [ start - 1]).upper()
        alt = str(genome["chr" + str(var["chrom"])]\
            [ start - 1: end]).upper()
        
        pos = str(start)
        
        return pos, ref, alt
    
    
    
    @staticmethod
    def del_genomic(var, genome):
        """
        Function
        --------
        Gets the genomic variant from the HGVS variant for del variants.
        
        Parameters
        ----------
        var : Variant dict.
        genome : The genome annotation Fasta.

        Returns
        -------
        pos : Genomic position of the mutation.
        ref : Ref at that genomic position.
        alt : Alt that has occured due to the mutation.

        """
        
        start = var["vep"]["start"] 
        end = var["vep"]["end"]
        ref = str(genome["chr" + str(var["chrom"])]\
            [ start - 2 : end ]).upper() # get the ref 
        alt = ref[0].upper()
        
        pos = str(start - 1)
        
        return pos, ref, alt
    
    
    
    @classmethod
    def snv_genomic(cls, var):
        """
        Function
        --------
        Gets the genomic variant from the HGVS variant for SNV variants.
        
        Parameters
        ----------
        cls : the class.
        var : variant dict.

        Returns
        -------
        pos : Genomic position of the mutation.
        ref : Ref at that genomic position.
        alt : Alt that has occured due to the mutation.

        """
        
        while True:    
            try:
                var_c = cls.hp.parse_hgvs_variant(var["hgvs"])
                var_g = cls.am.c_to_g(var_c) 
                
            except he.HGVSDataNotAvailableError:
                match = re.search("(^\w+\.\d+)(.*)", var["hgvs"])
                NM = str(match.group(1))
                newversion = int(NM[-1]) + 1
                NM = NM[0:-1] + str(newversion)[-1]
                var["hgvs"] = NM + str(match.group(2))
                var["warn"] = "HGVS error: genomic variant conversion was manually performed"
                continue
            break
        
        pos = str(var_g.posedit.pos.start.base)
        ref = str(var_g.posedit.edit.ref)
        alt = str(var_g.posedit.edit.alt)
        
        return pos, ref, alt
        
    
    
    ##########################################################################
    #                                   SPLICEAI                             #
    ##########################################################################     

    @staticmethod
    def get_spliceai(var_g, grch):
        
        # Define the url
        if "37" in grch:
            spliceaigrch = "hg=37"
        else:
            spliceaigrch = "hg=38"
            
        url = 'https://spliceailookup-api.broadinstitute.org/spliceai/?' +\
            spliceaigrch + '&variant=' + var_g
        
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
    
        # Make a dict of the scores
        scores = {
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
        
        return scores, mrnapos
    
    
    
    ##########################################################################
    #                             LOCATION METHODS                           #
    ##########################################################################     
    
    @classmethod
    def get_location(cls, transcripts, pos):
        """
        Function
        --------
        Gets the location of the variant (which intron or exon). It does so
        by getting the location of the variant of all the available ENSEMBL
        transcripts, and then selecting the most common one, to be as accurate
        as possible.
        
        Parameters
        ----------
        cls : the class.
        transcripts : list containing all of the ENSEMBL transcript IDs.
        pos : genomic position of the mutation.
        
        Raises
        ------
        RuntimeError : When no valid location is found.

        Returns
        -------
        location : Final selected location of the mutation.

        """
        
        location_counter = {}
        location_list = []
        for transcript in transcripts:
            location = cls.location_finder(transcript, pos)
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
            raise RuntimeError("Location not found for ")
        location = next(item for item in location_list if item["number"] == location_number)
        
        return location
    
    
        
    @staticmethod
    def location_finder(transcript, pos):
        """
        Function
        --------
        Finds the location of the mutation within a transcript. Does so by 
        comparing the genomic position with the exon start/ends. 

        Parameters
        ----------
        transcript : ENSEMBL transcript ID.
        pos : genomic position of the mutation.

        Returns
        -------
        location : Location of the position inside the transcript. Contains 
        which intron/exon and distance to either end.

        """
        
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
    #                           PREDICTION METHODS                           #
    ##########################################################################
    
    @classmethod
    def get_predictions(cls, var):
        """
        Function
        --------
        Controlling method to get the prediction for all SPLICEAI delta score 
        cut-offs. 

        Parameters
        ----------
        cls : The class.
        var : variant dict.

        Returns
        -------
        predictions : Dict of predictions for all cut-offs.

        """
        
        predictions = {
            0.8: "-",
            0.5: "-",
            0.2: "-",
            0.05: "-"
            }
        
        # Loop through the cutoffs
        maxcut = 1.0
        for cutoff in predictions:
            if sum(i >= cutoff and i < maxcut for i in var["scores"].values()):
                predictions[cutoff] = cls.predict(var, cutoff)
            maxcut = cutoff
            
        return predictions
    
    
    
    @classmethod
    def predict(cls, var, threshold):
        """
        Function
        --------
        Tries to predict a variants effect based on the SPLICEAI delta scores,
        the location of the variant and the position to start and end of the 
        location. 

        Parameters
        ----------
        var : variant dict.
        threshold : SPLICEAI delta score cut-off.

        Returns
        -------
        prediction : The final prediction of the variant/cut-off combination.

        """
            
        # Make it a bit simpler to work with
        acc_gain = var["scores"]["acceptor_gain"]
        acc_loss = var["scores"]["acceptor_loss"]
        don_gain = var["scores"]["donor_gain"]
        don_loss = var["scores"]["donor_loss"]
        
        
        
        # Acceptor gain
        if acc_gain >= threshold:
            
            # Acceptor gain alone
            if acc_loss < threshold and don_gain < threshold and don_loss < threshold:
                prediction = cls.gain(var, "acceptor")
                
            # Acceptor gain + loss
            if don_gain < threshold and don_loss < threshold:
                prediction = cls.gain_loss(var, "acceptor")
            
            # Acceptor gain + loss + donor loss/gain
            else:
                prediction = "complex"
                        
                
                
        # Acceptor loss
        elif acc_loss >= threshold:
            
            # Acceptor loss alone
            if don_gain <= threshold and don_loss <= threshold:
                
                prediction = cls.loss(var, "acceptor")

            # Acceptor loss + donor loss
            elif don_loss >= threshold and don_gain <= threshold:
                prediction = cls.loss_loss(var)

            # Acceptor loss + donor loss + donor gain
            else:
                prediction = "complex"
    
    
    
        # Donor gain
        elif don_gain >= threshold:
            
            # Donor gain alone
            if don_loss < threshold:
                prediction = cls.gain(var, "donor")
                
            # Donor gain + loss
            else:
                prediction = cls.gain_loss(var, "donor")
            
                
        # Donor loss
        elif don_loss >= threshold:
            
            # Donor loss alone
            prediction = cls.loss(var, "donor")
        
        return prediction
    
    
    
    @staticmethod
    def gain(var, site):
        """
        Function
        --------
        This method predicts gain-type mutations: donor or acceptor gain.
        Donor and acceptor gain mutation have a similar effect:either a small 
        insertion or deletion. Donor sites are at the 5' end of exons, and 
        acceptor sites at the 3' end. 
        
        Therefore, these mutations have the opposite effect of each other. 
        When an acceptor gain occurs upstream of the normal acceptor site, this
        will results in a deletion. But when a donor gain occurs upstream of 
        the normal donor site, an insertion occurs. 
        
        For that reason, upstream donor gain and downstream acceptor gain have
        the same effect on the transcript, as do downstream donor gain and 
        upstream acceptor gain.
        
        Parameters
        ----------
        var : class of the variant.
        site : site of the mutation; either donor or acceptor.

        Returns
        -------
        prediction : effect prediction of the mutation.

        """
        
        # Upstream
        if var["location"]["region"] == "exon" and site == "acceptor" or\
           var["location"]["region"] == "intron" and site == "donor":
            
            bp = var["mrnapos"]["acceptor_gain"] + var["location"]["start"] 
         
        # Downstream
        else:
            bp = var["mrnapos"]["acceptor_gain"] - var["location"]["start"] 
        
        # Upstream acceptor gain or downstream donor gain = del
        if site == "acceptor" and bp > 0 or \
           site == "donor" and bp < 0:
            effect = "del"
        
        # Downstream acceptor gain and upstream donor gain = ins
        else:
            effect = "ins"
        
        # Forgot why intronic acceptor sites are the next exon
        if var["location"]["region"] == "intron" and site == "acceptor":
            prediction = "Exon " + str(var["location"]["number"] + 1) + " " + \
                str(abs(bp)) + " bp " + effect
            
        else:
            prediction = "Exon " + str(var["location"]["number"]) + " " + \
                str(abs(bp)) + " bp " + effect
            
        return prediction
    
    
    
    @staticmethod
    def loss(var, site):
        """
        Function
        --------
        This method predict solo loss-type mutations: donor or acceptor loss.
        Donor or acceptor loss mutations lead to exon skips. 

        Parameters
        ----------
        var : class of the variant.
        site : site of the mutation; either donor or acceptor.

        Returns
        -------
        prediction : effect prediction of the mutation.

        """
        if var["location"]["region"] == "intron" and site == "acceptor":
            prediction = "Exon " + str(var["location"]["number"] + 1) + " skip"
            
        else:
            prediction = "Exon " + str(var["location"]["number"]) + " skip"
        
        return prediction
    
    
    
    @staticmethod
    def gain_loss(var, site):
        """
        Function
        --------
        This method predicts gain+loss-type mutations.
        Gain/loss-type mutations are very similar to gain mutations. They are
        still seperated in another method for clarity and because the basepairs
        are calculated differently.
        
        Parameters
        ----------
        var : class of the variant.
        site : site of the mutation; either donor or acceptor.

        Returns
        -------
        prediction : effect prediction of the mutation.

        """

        if site == "acceptor":
            bp = var["mrnapos"]["acceptor_gain"] - var["mrnapos"]["acceptor_loss"]
         

        else:
            bp = var["mrnapos"]["donor_gain"] - var["mrnapos"]["donor_loss"]
        
        # Upstream acceptor gain or downstream donor gain = del
        if site == "acceptor" and bp > 0 or \
           site == "donor" and bp < 0:
            effect = "del"
        
        # Downstream acceptor gain and upstream donor gain = ins
        else:
            effect = "ins"
        
        # Forgot why intronic acceptor sites are the next exon
        if var["location"]["region"] == "intron" and site == "acceptor":
            prediction = "Exon " + str(var["location"]["number"] + 1) + " " + \
                str(abs(bp)) + " bp " + effect
            
        else:
            prediction = "Exon " + str(var["location"]["number"]) + " " + \
                str(abs(bp)) + " bp " + effect
            
        return prediction
    
    
    
    @staticmethod
    def loss_loss(var):
        """
        Function
        --------
        This method predicts double loss mutations; albeit unlikely to occur.
        Double loss mutations can either lead to an exon skip, double exon skip
        or intron retention, based on the location of the mutations.
        

        Parameters
        ----------
        var : class of the variant.
        site : site of the mutation; either donor or acceptor.

        Returns
        -------
        prediction : effect prediction of the mutation.

        """

         # Exon skip if donor loss is downstream (positive 3') of acceptor loss
        if var["mrnapos"]["donor_loss"] > var["mrnapos"]["acceptor_loss"]:
            prediction = "Exon " + str(var["location"]["number"]) + " skip"
        
        # Double exon skip or intron retention if donor loss is upstream (negative 5') of acceptor loss
        elif var["mrnapos"]["donor_loss"] < var["mrnapos"]["acceptor_loss"]:
            prediction = "Exon " + str(var["location"]["number"] - 1) + \
                " & " + str(var["location"]["number"])  + " skip" + \
                        " OR intron " + str(var["location"]["number"] - 1) + \
                            " retention"
                                    
        # Exon skip if donor loss is downstream (positive 3')
        else:
            prediction = "Exon " + str(var["location"]["number"]) + " skip" 
            
        return prediction
    