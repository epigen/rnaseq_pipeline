## RSEQC


rule rseqc_gtf2bed:
    input:
        os.path.join(resource_path,"genome.gtf"),
    output:
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
        db=temp(os.path.join(result_path,"qc/rseqc/annotation.db")),
    log:
        "logs/rseqc_gtf2bed.log",
    conda:
        "../envs/gffutils.yaml"
    script:
        "../scripts/gtf2bed.py"


rule rseqc_junction_annotation:
    input:
        bam=os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.junctionanno.junction.bed"),
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_annotation/{sample}.log",
    params:
        extra=r"-q 255",  # STAR uses 255 as a score for unique mappers
        prefix=lambda w, output: output[0].replace(".junction.bed", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_annotation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log[0]} 2>&1"


rule rseqc_junction_saturation:
    input:
        bam=os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf"),
    priority: 1
    log:
        "logs/rseqc/rseqc_junction_saturation/{sample}.log",
    params:
        extra=r"-q 255",
        prefix=lambda w, output: output[0].replace(".junctionSaturation_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "junction_saturation.py {params.extra} -i {input.bam} -r {input.bed} -o {params.prefix} "
        "> {log} 2>&1"


rule rseqc_stat:
    input:
        os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.stats.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_stat/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "bam_stat.py -i {input} > {output} 2> {log}"

rule rseqc_infer:
    input:
        bam=os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.infer_experiment.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_infer/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "infer_experiment.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_innerdis:
    input:
        bam=os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_innerdis/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".inner_distance.txt", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "inner_distance.py -r {input.bed} -i {input.bam} -o {params.prefix} > {log} 2>&1"


rule rseqc_readdis:
    input:
        bam=os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
        bed=os.path.join(result_path,"qc/rseqc/annotation.bed"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.readdistribution.txt"),
    priority: 1
    log:
        "logs/rseqc/rseqc_readdis/{sample}.log",
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_distribution.py -r {input.bed} -i {input.bam} > {output} 2> {log}"


rule rseqc_readdup:
    input:
        os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.readdup.DupRate_plot.pdf"),
    priority: 1
    log:
        "logs/rseqc/rseqc_readdup/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".DupRate_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_duplication.py -i {input} -o {params.prefix} > {log} 2>&1"


rule rseqc_readgc:
    input:
        os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
    output:
        os.path.join(result_path,"qc/rseqc/{sample}.readgc.GC_plot.pdf"),
    priority: 1
    log:
        "logs/rseqc/rseqc_readgc/{sample}.log",
    params:
        prefix=lambda w, output: output[0].replace(".GC_plot.pdf", ""),
    conda:
        "../envs/rseqc.yaml"
    shell:
        "read_GC.py -i {input} -o {params.prefix} > {log} 2>&1"


rule multiqc:
    input:
        expand(
            os.path.join(result_path,"star/{sample}/Aligned.sortedByCoord.out.bam"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.junctionanno.junction.bed"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.junctionsat.junctionSaturation_plot.pdf"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.infer_experiment.txt"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.stats.txt"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.inner_distance_freq.inner_distance.txt"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.readdistribution.txt"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.readdup.DupRate_plot.pdf"),
            sample=list(samples.keys()),
        ),
        expand(
            os.path.join(result_path,"qc/rseqc/{sample}.readgc.GC_plot.pdf"),
            sample=list(samples.keys()),
        ),
        expand(
            "logs/rseqc/rseqc_junction_annotation/{sample}.log",
            sample=list(samples.keys()),
        ),
    output:
        os.path.join(result_path,"qc/multiqc_report.html"),
    log:
        "logs/multiqc.log",
    wrapper:
        "v3.5.3/bio/multiqc"
