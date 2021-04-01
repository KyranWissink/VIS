# Version 0.2 of VIpy

**Required file structure**
This version includes pyVEP
Some files are required for the script to work:

input
  |- hg37.fa or hg38.fa (genome fasta)
  |- variants.txt (txt of hgvs variants separated by newlines)
  |- preferred_transcripts.txt (OPTIONAL text file of preferred transcripts for NM numbers)
  
Script was tested using hg38; if hg37 does not work, that is probably why.
