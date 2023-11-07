## Containers instructions for qcflow-rnaseq

The pipeline is built suing anaconda.     

### Docker

Make sure you are runnning in the same directory as the Dockerfile.      

```
docker build . -t qcflow-rnaseq:v0.0.5
singularity build qcflow-rnaseq_v0.0.5.sif docker-daemon://qcflow-rnaseq:v0.0.5
```
This will create the singularity image to be used for the pipeline.     








