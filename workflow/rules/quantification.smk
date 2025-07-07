
# quantify number of reads per gene across samples
rule count_matrix:
    input:
        reads_per_gene = expand(
            os.path.join(result_path,"star/{sample}/ReadsPerGene.out.tab"),
            sample=list(samples.keys()),
        ),
    output:
       counts = os.path.join(result_path,"counts","counts.csv"),
    log:
        "logs/rules/count_matrix.log",
    params:
        samples=list(samples.keys()),
        strand=get_strandedness(annot_samples),
    resources:
        mem_mb=config.get("mem", "16000"),
    conda:
        "../envs/pandas.yaml"
    script:
        "../scripts/count_matrix.py"

# annotate genes using biomaRt & GTF
rule annotate_genes:
    input:
        counts = os.path.join(result_path,"counts","counts.csv"),
        gtf = os.path.join(resource_path,"genome.gtf"),
        fasta = os.path.join(resource_path,"genome.fasta"),
    output:
        gene_annotation=os.path.join(result_path,"counts","gene_annotation.csv"),
    params:
        species=get_bioc_species_name(),
        version=config["ref"]["release"],
    log:
        "logs/rules/annotate_genes.log",
    conda:
        "../envs/biomart.yaml"
    script:
        "../scripts/annotate_genes.R"

# create sample annotation file based on MultiQC general stats
rule sample_annotation:
    input:
        multiqc_stats = os.path.join(result_path, "report", "multiqc_report_data"),
    output:
        sample_annot = os.path.join(result_path, "counts", "sample_annotation.csv"),
    resources:
        mem_mb="4000",
    log:
        "logs/rules/sample_annotation.log"
    run:
        multiqc_df = pd.read_csv(os.path.join(input.multiqc_stats,"multiqc_general_stats.txt"), delimiter='\t', index_col=0).loc[list(samples.keys()),:]
        # merge by sample names (index) and drop redundant or unnecessary columns
        annot_df = pd.merge(annot_samples, multiqc_df, left_index=True, right_index=True, how='inner').drop(['bam_file', 'sample_name'], axis=1)
        # make column names R compatible
        annot_df.columns = (
                annot_df.columns
                  .str.replace(r'[^0-9A-Za-z_.]+', '_', regex=True) # replace any non-alphanumeric, non-underscore, non-dot with _
                  .str.replace(r'_+', '_', regex=True) # collapse multiple underscores
                  .str.strip('_') # strip leading/trailing underscores
                  .str.replace(r'^(\d)', r'X\1', regex=True) # if name starts with digit, prefix with 'X'
            )
        annot_df.to_csv(output.sample_annot)
