# Get all references and prepare them accordingly for downstream processing

rule get_genome:
    output:
        os.path.join(resource_path,"genome.fasta"),
    log:
        "logs/rules/get_genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
    cache: True
    wrapper:
        "v3.5.3/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        os.path.join(resource_path,"genome.gtf"),
    params:
        species=config["ref"]["species"],
        fmt="gtf",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        flavor="",
    cache: True
    log:
        "logs/rules/get_annotation.log",
    wrapper:
        "v3.5.3/bio/reference/ensembl-annotation"


rule genome_faidx:
    input:
        os.path.join(resource_path,"genome.fasta"),
    output:
        os.path.join(resource_path,"genome.fasta.fai"),
    log:
        "logs/rules/genome_faidx.log",
    cache: True
    wrapper:
        "v3.5.3/bio/samtools/faidx"


rule bwa_index:
    input:
        os.path.join(resource_path,"genome.fasta"),
    output:
        multiext(os.path.join(resource_path,"genome.fasta"), ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/rules/bwa_index.log",
    resources:
        mem_mb=369000,
    cache: True
    wrapper:
        "v3.5.3/bio/bwa/index"


rule star_index:
    input:
        fasta=os.path.join(resource_path,"genome.fasta"),
        annotation=os.path.join(resource_path,"genome.gtf"),
    output:
        directory(os.path.join(resource_path,"star_genome")),
    resources:
        mem_mb=4*int(config.get("mem", "16000")),
    threads: 8
    params:
        extra=lambda wc, input: f"--sjdbGTFfile {input.annotation} --sjdbOverhang 100",
    log:
        "logs/rules/star_index_genome.log",
    cache: True
    wrapper:
        "v3.5.3/bio/star/index"
