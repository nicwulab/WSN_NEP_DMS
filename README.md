# NEP Deep mutation scanning and IAV RNA species(mRNA/cRNA/vRNA) quantification

## Dependencies ##
* python=3.6
* [snakemake](https://snakemake.readthedocs.io/en/stable/)
* [flash](https://github.com/dstreett/FLASH2)
* [seqtk](https://github.com/lh3/seqtk)
* [cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)
* [STAR](https://github.com/alexdobin/STAR)
* [pysam](https://github.com/pysam-developers/pysam)
## Installation ##
Install everything dependencies by conda:

```conda create -n NEP -c bioconda -c anaconda python=3.6 seqtk flash```

Before running the analysis, do:

```
source activate NEP(for Mac)
```
## DMS analysis ##
1. set the env
2. Run analysis:
    - go into ```/script/``` folder
    - set the ```PROJECT_PATH``` variable in ```nep_pipeline.smk``` file accordingly
    - do: ```snakemake --use-conda --cores 4 -s nep_pipeline.smk -j 2``` to excute the analysis
   

## IAV RNA species analysis ##


### 1.1 make reference
```
RNAseq_script/make_ref.sh
```
### 1.2 align to IAV genome by STAR
```
RNAseq_script/pipeline_align.sh
```
### 1.3 RNA quantification (htseq-count)
```
RNAseq_script/pipeline_count.smk
python RNAseq_script/align2count.py
```
### 1.4 Deletion analysis
```
python RNAseq_script/align2DIs.py
```
### 1.5 cRNA and mRNA ratio estimation by 3 end
```
python RNAseq_script/ratio_by3prime_v2.py
```