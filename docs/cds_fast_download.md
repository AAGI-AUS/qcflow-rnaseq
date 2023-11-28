# Download fasta files

## Download singularity image
Link to [image](https://hub.docker.com/r/staphb/ncbi-datasets/tags)

```
singularity pull docker://staphb/ncbi-datasets:latest
```

## Download data

Keep the fasta files separated whether you want to test for different types of contamination.      

```
TAXONID=4890
singularity run ncbi-datasets_latest.sif datasets download genome taxon $TAXONID  --include cds

unzip ncbi_dataset.zip

find -name "*fna" -exec cat {} >> taxid${TAXONID}.fasta \;
```

Useful taxon ids 

|Used | Taxon  | Name            |
|-|--------|-----------------|
|*|4890   | Ascomycota      |
| |147541 | Dothideomycetes |
| |55067  | Phaeosphaeria   |
| |4565   | wheat           |
|*|4513   | barley          |
|*|6946   | Mites           |
|*|30262  | Thrips          |
|*|27482  | Aphids          |

## Additional contamination screenings

* Plastid genomes download from [here](https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/)
* UniVec [database](https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/)
