## VIpy: HGVS variant interpretation using SPLICEAI<br>
This script allows direct HGVS mutation variant prediction using SpliceAI. <br>
<br>
### Prerequisites
This script requires genome annotation for the genome the user provides. These can be downloaded here:<br>
hg38: http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz<br>
hg19: http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz<br>
<br>
This script requires some dependencies to run. These can be found on their respective GitHub pages:<br>
https://github.com/biocommons/hgvs<br>
https://github.com/mdshw5/pyfaidx<br>
https://github.com/gawbul/pyEnsemblRest<br>
https://github.com/pandas-dev/pandas<br>
<br>
Alternatively, these can be installed directly via:<br>
```sh
pip install hgvs
pip install pyfaidx
pip install pyensemblrest
pip install pandas
```

### Usage
The script can be run directly from the command line:<br>
```sh
HGVSpredict.py -I input -O output -G genome -P preferred_transcript (optional)
```
<br>
<br>
![RUG logo](https://www.rug.nl/about-ug/practical-matters/huisstijl/huisstijl-basiselementen/images/rugr_logonl_rood_rgb-web.png)<br><br>
**Kyran Wissink**<br>Student Biomedical Sciences<br>University of Groningen<br>github.com/KyranWissink<br>k.wissink@student.rug.nl
