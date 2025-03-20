
# quantify number of reads per gene across samples
rule count_matrix:
    input:
        expand(
            os.path.join(result_path,"star/{sample}/ReadsPerGene.out.tab"),
            sample=list(samples.keys()),
        ),
    output:
       os.path.join(result_path,"counts/counts.csv"),
    log:
        "logs/count_matrix.log",
    params:
        samples=list(samples.keys()),
        strand=get_strandedness(annot_samples),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count_matrix.py"


# currently unused
# rule gene_2_symbol:
#     input:
#         counts="{prefix}.csv",
#     output:
#         symbol="{prefix}.symbol.csv",
#     params:
#         species=get_bioc_species_name(),
#     log:
#         "logs/gene2symbol/{prefix}.log",
#     conda:
#         "../envs/biomart.yaml"
#     script:
#         "../scripts/gene2symbol.R"
