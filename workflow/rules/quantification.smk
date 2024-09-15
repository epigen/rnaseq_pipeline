rule count_matrix:
    input:
        expand(
            "results/star/{unit.sample_name}_{unit.unit_name}/ReadsPerGene.out.tab",
            unit=units.itertuples(),
        ),
    output:
        "results/counts/all.tsv",
    log:
        "logs/count-matrix.log",
    params:
        samples=units["sample_name"].tolist(),
        strand=get_strandedness(units),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count-matrix.py"


rule gene_2_symbol:
    input:
        counts="{prefix}.tsv",
    output:
        symbol="{prefix}.symbol.tsv",
    params:
        species=get_bioc_species_name(),
    log:
        "logs/gene2symbol/{prefix}.log",
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/gene2symbol.R"
