# NEP Deep mutation scanning

## Dependencies ##
* python=3.6
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
## Installation ##
Install everything dependencies by conda:

```conda create -n NEP -c bioconda -c anaconda python=3.6 seqtk flash```

Before running the analysis, do:

```
source activate SARS(for Mac)
```
## Steps ##
1. set the env
2. Run analysis:
    - go into ```/script/``` folder
    - set the ```PROJECT_PATH``` variable in ```nep_pipeline.smk``` file accordingly
    - do: ```snakemake -s pipeline.smk``` to excute the analysis