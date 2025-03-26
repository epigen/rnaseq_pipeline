
# libraries
import sys
import pandas as pd

# logging
sys.stderr = open(snakemake.log[0], "w")

#### config

# input
reads_per_gene_list = snakemake.input

# output
counts_path = snakemake.output["counts"]

# params
samples = snakemake.params.samples
strand = snakemake.params.strand

# get column based on strandedness
def get_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 1  # non stranded protocol
    elif strandedness == "yes":
        return 2  # 3rd column
    elif strandedness == "reverse":
        return 3  # 4th column, usually for Illumina truseq
    else:
        raise ValueError(
            (
                "'strandedness' column should be empty or have the "
                "value 'none', 'yes' or 'reverse', instead has the "
                "value {}"
            ).format(repr(strandedness))
        )

# aggregate and save counts across samples (read_per_genes files)
counts = [
    pd.read_table(
        f, index_col=0, usecols=[0, get_column(strandedness)], header=None, skiprows=4
    )
    for f, strandedness in zip(reads_per_gene_list, strand)
]

for t, sample in zip(counts, samples):
    t.columns = [sample]

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
matrix.to_csv(counts_path)
