# Version 0.3 of VIpy<Br>

**Command line arguments**
arg1: variants in txt file
arg2: hg19 OR hg38

This version does **NOT** requires pyVEP
**Required file structure**<br>
Some files are required for the script to work:<br>

input (folder name)<br>
  |- hg37.fa or hg38.fa (genome fasta)<br>
  |- variants.txt (txt of hgvs variants separated by newlines)<br>
  |- preferred_transcripts.txt (OPTIONAL text file of preferred transcripts for NM numbers)<br><br>
  
Script was tested using hg38 and hg37
