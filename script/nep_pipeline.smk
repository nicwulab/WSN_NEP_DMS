# only variable needed to change
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import numpy as np

PROJECT_PATH='/Users/yiquan/PycharmProjects/NEP Project'
R1 = PROJECT_PATH + '/data/{SAMPLENAME}_L001_R1_001.fastq.gz'
R2 = PROJECT_PATH + '/data/{SAMPLENAME}_L001_R2_001.fastq.gz'
SAMPLENAMES, = glob_wildcards(R1)


rule all:
    input:
        expand('{sample}.tsv', sample=SAMPLES)

rule merge_RMTS:
    input:
        FQ1=R1,
        FQ2=R2
    output:
        COMBINED_FAS
    shell:
        'python Fastq2ErrorFreeFasta.py -i {input} -o {output} -b 0-0 -p 1-8 -d 2 -F _R1_ -R _R2_ -e 0.9 -s 2'

rule fas2fq:
    input:
        COMBINED_FAS
    output:
        COMBINED_FQ
    shell:
        'seqtk seq -F '#' {input} > {output}'

rule flash:
    input:
        COMBINED_FQ
    output:
        ASSEMBLED_FQ
    shell:
        'flash -m70 -M 80 --interleaved-input {input} -o {output}'

rule fq2count:
    input:
        ASSEMBLED_FQ
    output:
        TABLE
    shell:
        'python fastq2count.py {input} {output}'


