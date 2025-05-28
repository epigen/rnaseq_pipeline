
# check for each sample's bam files if the provided read type (single or paired) is correct
rule check_read_type:
    input:
        bams = lambda wc: annot.loc[wc.sample, "bam_file"],
    output:
        check = os.path.join(result_path,".check_read_type","{sample}.done"),
    params:
        read_type = lambda wc: 'SE' if samples[wc.sample]['read_type'] == 'single' else 'PE',
        samtools_threads = lambda wc, threads: int(threads) - 1,
    threads: 10
    resources:
            mem_mb=config.get("mem", "16000"),
    log:
        "logs/rules/check_read_type_{sample}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        for bam_file in {input.bams}; do
            # Use samtools to count reads flagged as paired (flag 0x1)
            paired_count=$(samtools view --threads {params.samtools_threads} -c -f 0x1 "$bam_file")

            # If any reads have the paired flag, we consider the file paired-end.
            if [ "$paired_count" -gt 0 ]; then
              actual_type="PE"
            else
              actual_type="SE"
            fi

            # Compare the detected type with the expected type
            if [ "$actual_type" != "{params.read_type}" ]; then
              echo "Error: BAM file type ($actual_type) does not match expected type {params.read_type}. Exiting."
              exit 1
            else
              echo "BAM file type matches expected type."
            fi
        done

        touch {output.check}
        """

# Merge uBAM files, convert to interleaved FASTQ, trim and filter using fastp, then de-interleave for alignment
# de-interleaving from here: https://gist.github.com/nathanhaigh/3521724
# tested by comparing the output to seqfu interleave (https://telatin.github.io/seqfu2/tools/deinterleave.html)
rule trim_filter:
    input:
        bams = lambda wc: annot.loc[wc.sample, "bam_file"],
        read_type_check = os.path.join(result_path,".check_read_type","{sample}.done"),
        adapter_fasta = config["adapter_fasta"] if config["adapter_fasta"]!="" else [],
    output:
        fastq_filtered_R1 = temp(os.path.join(result_path,"fastp","{sample}","{sample}_R1.filtered.fastq.gz")),
        fastq_filtered_R2 = temp(os.path.join(result_path,"fastp","{sample}","{sample}_R2.filtered.fastq.gz")),
        fastp_html = os.path.join(result_path,"fastp","{sample}","{sample}.fastp.html"),
        fastp_json = os.path.join(result_path,"fastp","{sample}","{sample}.fastp.json"),
    params:
        read_type = lambda wc: 'SE' if samples[wc.sample]['read_type'] == 'single' else 'PE',
        # samtools fastq args
        fastq_opts = lambda wc: "-N" if samples[wc.sample]['read_type'] == 'paired' else "",
        samtools_threads = lambda wc, threads: int(threads) - 1,
        # fastp adapter trimming and filtering args
        fastp_args = config["fastp_args"] if config["fastp_args"] != "" else "",
        adapter_fasta = "--adapter_fasta " + config["adapter_fasta"] if config["adapter_fasta"] !="" else "",
        interleaved_in = lambda wc: "--interleaved_in" if samples[wc.sample]['read_type'] == 'paired' else "",
    threads: 10
    resources:
            mem_mb=config.get("mem", "16000"),
    log:
        samtools = "logs/samtools/{sample}.log",
        fastp = "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yaml"
    shell:
        """
        samtools merge --threads {params.samtools_threads} -u - {input.bams} 2>> "{log.samtools}" | \
        samtools fastq --threads {params.samtools_threads} {params.fastq_opts} - 2>> "{log.samtools}" | \
        fastp {params.fastp_args} {params.adapter_fasta} {params.interleaved_in} --thread {threads} --stdin --stdout  --html "{output.fastp_html}" --json "{output.fastp_json}" 2> "{log.fastp}" | \
        {{
          if [ "{params.read_type}" = "PE" ]; then
              # For paired-end: de-interleave the FASTQ output and compress R1 and R2
              paste - - - - - - - - | tee >(cut -f 1-4 | tr "\\t" "\\n" | pigz --best --processes {threads} > "{output.fastq_filtered_R1}") | cut -f 5-8 | tr "\\t" "\\n" | pigz --best --processes {threads} > "{output.fastq_filtered_R2}"
          else
              # For single-end: compress output for R1 and create an empty dummy file for R2
              pigz --best --processes {threads} > "{output.fastq_filtered_R1}"
              touch "{output.fastq_filtered_R2}"
          fi
        }}
        """

# align reads directly from temporary, trimmed and filtered gzipped FASTQ files
rule align:
    input:
        fastq_filtered_R1 = os.path.join(result_path,"fastp","{sample}","{sample}_R1.filtered.fastq.gz"),
        fastq_filtered_R2 = os.path.join(result_path,"fastp","{sample}","{sample}_R2.filtered.fastq.gz"),
        index = os.path.join(resource_path,"star_genome"),
        gtf = os.path.join(resource_path,"genome.gtf"),
    output:
        bam = os.path.join(result_path,"star","{sample}","Aligned.sortedByCoord.out.bam"),
        bai = os.path.join(result_path,"star","{sample}","Aligned.sortedByCoord.out.bam.bai"),
        reads_per_gene = os.path.join(result_path,"star","{sample}","ReadsPerGene.out.tab"),
    resources:
        # dynamic memory allocation based on attempts (multiple attempts can be configured with --retries X)
        mem_mb=lambda wildcards, attempt: attempt*int(config.get("mem", "32000")),
    threads: 24
    log:
        "logs/star/{sample}.log",
    conda:
        "../envs/star.yaml"
    params:
        star_input = lambda wc, input: f'"{input.fastq_filtered_R1}"' if samples[wc.sample]['read_type'] == 'single' else f'"{input.fastq_filtered_R1}" "{input.fastq_filtered_R2}"',
        star_args = config['star_args'],
        result_dir = lambda wc: os.path.join(result_path,"star",f"{wc.sample}"),
        samtools_threads = lambda wc, threads: int(threads) - 1,
    shell:
        """
        # run STAR alignment
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

        # index BAM file
        samtools index --threads {params.samtools_threads} "{output.bam}" "{output.bai}"
        """
