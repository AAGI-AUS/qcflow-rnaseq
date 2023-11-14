## Info about Pawsey and running large jobs

### Genome indexing

* STAR genome indexing uses about 30GB for the human genome and about 44 GB for the barley genome (4.3 Gbp). Make sure you allocate the right memory with ```--mem```.     
* HISAT2 requires a large memory for indexing (especially using ```--ss``` and ```--exon```). Make sure you allocate enough physical memory. HISAT2 doesn't work with virtual memory. In this case use ```--mem-per-cpu``` option. As for large genomes make sure you run on the ```highmem``` partition of Pawsey.      

HISAT index generation can be used by modifying the indexing parameters such as ```--large-index```. If that doesn't work you can manually change parameters such as ```--noauto```, increasing to 8, 16, 32.. the ```--bmaxdivn``` parameter and changing ```--dcv``` to 4096. I have generated a function that automatically takes care of that.           
More details about the memory requirements can be found [here](https://www.biostars.org/p/9521274/).     
