# only variable needed to change


PROJECT_PATH='/Users/yiquan/PycharmProjects/NEP_Project'
FQ = PROJECT_PATH + '/fastq/{SAMPLENAME}/L001_R1_001.fastq.gz'
SAMPLENAMES, = glob_wildcards(FQ)
FQ_FOLDER = PROJECT_PATH + '/fastq/{SAMPLENAME}'
#print(SAMPLENAMES)
RESULT_PATH = PROJECT_PATH + '/results/{SAMPLENAME}'
RM_TAG_FAS = RESULT_PATH + '/rm_tag.fa'
RM_PRIMER_FAS = RESULT_PATH + '/rm_primer.fa'
UNTRIM_FAS = RESULT_PATH + '/untrim.fa'
RM_PRIMER_FQ = RESULT_PATH + '/rm_tag.fq'
ASSEMBLED_FQ = RESULT_PATH + '/assembled.fq'
TABLE = RESULT_PATH + '/nep_mut.tsv'

rule all:
    input:
        expand(TABLE, SAMPLENAME=SAMPLENAMES),
        expand(RM_TAG_FAS, SAMPLENAME = SAMPLENAMES),
        expand(RM_PRIMER_FAS, SAMPLENAME = SAMPLENAMES),
        expand(RM_PRIMER_FQ, SAMPLENAME = SAMPLENAMES),
        expand(ASSEMBLED_FQ, SAMPLENAME = SAMPLENAMES)

rule merge_RMTS:
    input:
        FQ_FOLDER
    output:
        RM_TAG_FAS
    conda:
        PROJECT_PATH+ "/env/merge_RMTS.yaml"
    shell:
        'python Fastq2ErrorFreeFasta.py -i {input} -o {output} -b 0-0 -p 1-8 -d 3 -F _R1_ -R _R2_ -e 0.7 -s 2'

# remove primer(only for paired-primer are removed) &filter out small reads(-m 100:100(R1:R2))
rule rm_primer:
    input:
        RM_TAG_FAS
    params:
        rename=lambda wc: "'{id}_{adapter_name}'"
    output:
        trim_o=RM_PRIMER_FAS,
        untrim_o = UNTRIM_FAS
    shell:
        '''cutadapt --pair-adapters --interleaved '''\
        '''-g "trim=GACGTCGAGGAGAATCCCGGGCCC;e=0.2"  '''\
        '''-G "trim=AAACAAGGGTGTTTTTTATTATTA;e=0.2" '''\
        '''--rename={params} '''\
        '''--untrimmed-output {output.untrim_o} '''
        '''-O 10 -m 100:100 -o {output.trim_o} {input}'''

rule fas2fq:
    input:
        RM_PRIMER_FAS
    output:
        RM_PRIMER_FQ
    shell:
        "seqtk seq -F '#' {input} > {output}"

rule flash:
    input:
        RM_PRIMER_FQ
    output:
        ASSEMBLED_FQ
    shell:
        "flash -m 70 -M 80 --interleaved-input {input} -c > {output} "

rule fq2count:
    input:
        ASSEMBLED_FQ
    output:
        TABLE
    shell:
        'python fastq2count.py {input} {output}'


