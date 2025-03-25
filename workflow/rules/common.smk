
# used in rule count_matrix
def get_strandedness(annot_samples):
    if "strandedness" in annot_samples.columns:
        return annot_samples["strandedness"].tolist()
    else:
        strand_list = ["none"]
        return strand_list * annot_samples.shape[0]

# used for rule annotate_genes
def get_bioc_species_name():
    first_letter = config["ref"]["species"][0]
    subspecies = config["ref"]["species"].split("_")[1]
    return first_letter + subspecies

