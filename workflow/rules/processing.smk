
# Merge uBAM files, convert to interleaved FASTQ, trim and filter using fastp, then de-interleave for alignment
# de-interleaving from here: https://gist.github.com/nathanhaigh/3521724
# tested by comparing to output from seqfu interleave (https://telatin.github.io/seqfu2/tools/deinterleave.html)
rule trim_filter:
    input:
        bams = lambda wc: annot.loc[wc.sample, "bam_file"],
    output:
        fastq_filtered_R1 = temp(os.path.join(result_path,"fastp","{sample}","{sample}_R1.filtered.fastq.gz")),
        fastq_filtered_R2 = temp(os.path.join(result_path,"fastp","{sample}","{sample}_R2.filtered.fastq.gz")),
        fastp_html = os.path.join(result_path,"fastp","{sample}","{sample}.fastp.html"),
        fastp_json = os.path.join(result_path,"fastp","{sample}","{sample}.fastp.json"),
        fastp_log = os.path.join(result_path,"fastp","{sample}","{sample}.fastp.log"),
        samtools_log = os.path.join(result_path,"fastp","{sample}","{sample}.samtools.log"),
    params:
        read_type = lambda wc: 'SE' if samples[wc.sample]['read_type'] == 'single' else 'PE',
        # samtools fastq args
        fastq_opts = lambda wc: "-N" if samples[wc.sample]['read_type'] == 'paired' else "",
        # fastp adapter trimming and filtering args
        interleaved_in = lambda wc: "--interleaved_in" if samples[wc.sample]['read_type'] == 'paired' else "",
        adapter_sequence = "-a " + config["adapter_sequence"] if config["adapter_sequence"] != "" else "",
        adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else "",
        fastp_args = config["fastp_args"] if config["fastp_args"] != "" else "",
    threads: 10
    resources:
            mem_mb=config.get("mem", "16000"),
    log:
        "logs/rules/trim_filter_{sample}.log"
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        samtools merge -u - "{input.bams}" 2>> "{output.samtools_log}" | \
        samtools fastq {params.fastq_opts} - 2>> "{output.samtools_log}" | \
        fastp {params.adapter_sequence} {params.adapter_fasta} {params.interleaved_in} --stdin --stdout {params.fastp_args} --html "{output.fastp_html}" --json "{output.fastp_json}" 2> "{output.fastp_log}" | \
        {{
          if [ "{params.read_type}" = "PE" ]; then
              # deinterleave and gzip into R1 and R2 fastq.gz
              paste - - - - - - - - | tee >(cut -f 1-4 | tr "\\t" "\\n" | pigz --best --processes {threads} > "{output.fastq_filtered_R1}") | cut -f 5-8 | tr "\\t" "\\n" | pigz --best --processes {threads} > "{output.fastq_filtered_R2}"
          else
              # gzip R1
              pigz --best --processes {threads} > "{output.fastq_filtered_R1}"
              # touch R2
              touch "{output.fastq_filtered_R2}"
          fi
        }}
        """

# align reads directly from trimmed and filtered BAM files
rule align:
    input:
        fastq_filtered_R1 = os.path.join(result_path,"fastp","{sample}","{sample}_R1.filtered.fastq.gz"),
        fastq_filtered_R2 = os.path.join(result_path,"fastp","{sample}","{sample}_R2.filtered.fastq.gz"),
        index = os.path.join(resource_path,"star_genome"),
        gtf = os.path.join(resource_path,"genome.gtf"),
    output:
        bam = os.path.join(result_path,"star","{sample}","Aligned.sortedByCoord.out.bam"),
        reads_per_gene = os.path.join(result_path,"star","{sample}","ReadsPerGene.out.tab"),
    resources:
        mem_mb=config.get("mem", "16000"),
    threads: 24
    log:
        "logs/star/{sample}.log",
    conda:
        "../envs/star.yaml"
    params:
        star_input = lambda wc, input: f'"{input.fastq_filtered_R1}"' if samples[wc.sample]['read_type'] == 'single' else f'"{input.fastq_filtered_R1}" "{input.fastq_filtered_R2}"',
        star_args = config['star_args'],
        result_dir = lambda wc: os.path.join(result_path,"star",f"{wc.sample}"),
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir "{input.index}" \
             --readFilesType Fastx \
             --readFilesCommand zcat \
             --readFilesIn {params.star_input} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode GeneCounts \
             --sjdbGTFfile "{input.gtf}" \
             {params.star_args} \
             --outFileNamePrefix {params.result_dir}/ \
             > {log} 2>&1
        """
