# only variable needed to change

PROJECT_PATH = '/home/yiquan2/NEP'
REF = PROJECT_PATH + '/data/ref'
ANNOT = REF + '/ANNOT.gff3'

# Use glob_wildcards to get sample names dynamically
SAMPLENAMES, = glob_wildcards(PROJECT_PATH + '/data/raw/{SAMPLENAME}_L00M_R1_001.fastq.gz')

RESULT_PATH = PROJECT_PATH + '/result/{SAMPLENAME}'
BAM = RESULT_PATH + '/{SAMPLENAME}Aligned.sortedByCoord.out.bam'
SORTED_BAM = RESULT_PATH + '/{SAMPLENAME}_sortedbyname.bam'
COUNTS = RESULT_PATH + '/{SAMPLENAME}_counts.txt'
BAI = RESULT_PATH + '/{SAMPLENAME}Aligned.sortedByCoord.out.bam.bai'  


rule all:
    input:
        expand(BAI, SAMPLENAME = SAMPLENAMES),
        expand(COUNTS, SAMPLENAME = SAMPLENAMES)

rule sort_bam:
    input:
        BAM=BAM
    output:
        SORTED_BAM=SORTED_BAM,
        BAI=BAI
    shell:
        "samtools sort -n -o {output.SORTED_BAM} {input.BAM}"
        " && samtools index {input.BAM}"

rule htseq_count:
    input:
        SORTED_BAM=SORTED_BAM
    output:
        COUNTS=COUNTS
    params:
        annot_file = ANNOT
    shell:
        "htseq-count {input.SORTED_BAM} {params.annot_file} > {output.COUNTS}"
