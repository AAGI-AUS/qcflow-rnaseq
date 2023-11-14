## Troubleshooting

### Tools
* ```rnaseqc``` and input ```*gff``` file: make sure the input file contains ```gene``` entries. ``gtf2bed.py``` allows to convert the gtf file to bed input.    

### System
* I receive the error ```error from nextflow (ERR): mkfifo(/tmp/17.inpipe1) failed.``` when using the latest version of nextflow ```nextflow/23.04.3```. Changing to an earlier version removes the error.    
