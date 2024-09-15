# rule align:
#     input:
#         unpack(get_fq),
#         index="resources/star_genome",
#         gtf="resources/genome.gtf",
#     output:
#         aln="results/star/{sample}_{unit}/Aligned.sortedByCoord.out.bam",
#         reads_per_gene="results/star/{sample}_{unit}/ReadsPerGene.out.tab",
#     log:
#         "logs/star/{sample}_{unit}.log",
#     params:
#         idx=lambda wc, input: input.index,
#         extra=lambda wc, input: f'--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile {input.gtf} {config["params"]["star"]}',
#     threads: 24
#     wrapper:
#         "v3.5.3/bio/star/align"

# align reads directly from mutliple BAM files
rule align:
    input:
        bams=get_raw_bams,
        index=os.path.join(resource_path,"star_genome"),
        gtf=os.path.join(resource_path,"genome.gtf"),
    output:
        bam=os.path.join(result_path,"star","{sample}","Aligned.sortedByCoord.out.bam"),
        reads_per_gene=os.path.join(result_path,"star","{sample}","ReadsPerGene.out.tab"),
    log:
        "logs/star/{sample}.log",
    conda:
        "../envs/star.yaml"
    params:
        bam_files=lambda wc, input: ','.join(input.bams),
        read_type=lambda wc: 'SE' if samples[wc.sample]['read_type'] == 'single' else 'PE',
        extra=config['params']['star'],
    threads: 24
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.index} \
             --readFilesType SAM {params.read_type} \
             --readFilesCommand samtools view -h \
             --readFilesIn {params.bam_files} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile {input.gtf} \
             {params.extra} \
             --outFileNamePrefix $(dirname {output.bam})/ \
             > {log} 2>&1
        """


