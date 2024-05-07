# NEP Deep mutation scanning

## Dependencies ##
* python=3.6
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
## Installation ##
Install everything dependencies by conda:

```conda create -n NEP -c bioconda -c anaconda python=3.6 seqtk flash```

Before running the analysis, do:

```
source activate NEP(for Mac)
```
## Steps ##
1. set the env
2. Run analysis:
    - go into ```/script/``` folder
    - set the ```PROJECT_PATH``` variable in ```nep_pipeline.smk``` file accordingly
    - do: ```snakemake --use-conda --cores 4 -s nep_pipeline.smk -j 2``` to excute the analysis
   
3. The excuted workflow is as following: 

![workflow](https://github.com/Wangyiquan95/NEP/blob/main/script/dag.pdf?raw=true)
