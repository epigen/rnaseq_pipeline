
def get_raw_bams(wildcards):
    return annot.loc[wildcards.sample, "bam_file"]

# used in rule count_matrix
def get_strandedness(annot_samples):
    if "strandedness" in annot_samples.columns:
        return annot_samples["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * annot_samples.shape[0]

def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if pd.isna(unit["fq1"]):
        # SRA sample (always paired-end for now)
        accession = unit["sra"]
        return expand("sra/{accession}_{read}.fastq", accession=accession, read=[1, 2])

    if unit["fq1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["fq2"]):
        # single end local sample
        return "pipe/cutadapt/{S}/{U}.fq1.fastq{E}".format(
            S=unit.sample_name, U=unit.unit_name, E=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/cutadapt/{S}/{U}.{{read}}.fastq{E}".format(
                S=unit.sample_name, U=unit.unit_name, E=ending
            ),
            read=["fq1", "fq2"],
        )


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob.glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0
    return files

# # used in get_fq() below
# def is_paired_end(sample):
#     sample_units = units.loc[sample]
#     fq2_null = sample_units["fq2"].isnull()
#     sra_null = sample_units["sra"].isnull()
#     paired = ~fq2_null | ~sra_null
#     all_paired = paired.all()
#     all_single = (~paired).all()
#     assert (
#         all_single or all_paired
#     ), "invalid units for sample {}, must be all paired end or all single end".format(
#         sample
#     )
#     return all_paired

# # used before in rule align
# def get_fq(wildcards):
#     if config["trimming"]["activate"]:
#         # activated trimming, use trimmed data
#         if is_paired_end(wildcards.sample):
#             # paired-end sample
#             return dict(
#                 zip(
#                     ["fq1", "fq2"],
#                     expand(
#                         "results/trimmed/{sample}_{unit}_{group}.fastq.gz",
#                         group=["R1", "R2"],
#                         **wildcards,
#                     ),
#                 )
#             )
#         # single end sample
#         return {
#             "fq1": "results/trimmed/{sample}_{unit}_single.fastq.gz".format(**wildcards)
#         }
#     else:
#         # no trimming, use raw reads
#         u = units.loc[(wildcards.sample, wildcards.unit)]
#         if pd.isna(u["fq1"]):
#             # SRA sample (always paired-end for now)
#             accession = u["sra"]
#             return dict(
#                 zip(
#                     ["fq1", "fq2"],
#                     expand(
#                         "sra/{accession}_{group}.fastq",
#                         accession=accession,
#                         group=["R1", "R2"],
#                     ),
#                 )
#             )
#         if not is_paired_end(wildcards.sample):
#             return {"fq1": f"{u.fq1}"}
#         else:
#             return {"fq1": f"{u.fq1}", "fq2": f"{u.fq2}"}



# used for rule gene_2_symbol
def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies

