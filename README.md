## VIpy: HGVS variant interpretation using SPLICEAI
This script allows direct HGVS mutation variant prediction using SpliceAI. 
<br>
This entire script is based on SpliceAI. The code can be found on their GitHub:<br>
https://github.com/Illumina/SpliceAI
<br>

### Prerequisites
This script requires genome annotation for the genome the user provides. These can be downloaded here:<br>
hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
<br>
<br>
This script requires some dependencies to run. These can be found on their respective GitHub pages:<br>
https://github.com/biocommons/hgvs
https://github.com/mdshw5/pyfaidx
https://github.com/gawbul/pyEnsemblRest
https://github.com/pandas-dev/pandas

Alternatively, these can be installed directly via:
```sh
pip install hgvs
pip install pyfaidx
pip install pyensemblrest
pip install pandas
```

### Usage
The script can be run directly from the command line:
```sh
HGVSpredict.py -I input -O output -G genome -P preferred_transcript (optional)
```

### Code flow
* Check arguments
* Validate variants with preferred transcripts
* Per-variant runs:
* * Conversion from HGVS to genomic variant 
* * Locating the mutation within the gene
* * Get SpliceAI scores
* * Predict transcript effect based on location and scores

![RUG logo](https://www.rug.nl/about-ug/practical-matters/huisstijl/huisstijl-basiselementen/images/rugr_logonl_rood_rgb-web.png)<br>
Kyran Wissink
Student Biomedical Sciences
University of Groningen
github.com/KyranWissink
k.wissink@student.rug.nl
