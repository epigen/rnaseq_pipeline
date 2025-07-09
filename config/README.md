# Configuration

You need one configuration file (`config.yaml`) and one annotation file (`annotation.csv`) to run the complete workflow. You can use the provided examples as starting point. If in doubt read the comments in the config and/or try the default values.

- project configuration (`config/config.yaml`): different for every project/dataset and configures the processing and quantification. The fields are described within the file.
- annotation file (`annotation.csv`): `CSV` file consisting of one technical sequencing unit per row (i.e., one sample can include multiple sequencing units e.g., if you have several runs or lanes per sample, hence mutliple rows per sample) and 4 mandatory columns:
  - sample_name: No whitespace, special characters, starting with numbers or hyphen (`-`) etc. allowed. We recommend `snake_case`.
  - read_type: "single" or "paired"
  - bam_file: path to the raw/unaligned/unmapped [uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format) files. No whitespaces in file paths of names allowed.
  - strandedness: To get the correct `geneCounts` from `STAR` output, you can provide information on the strandedness of the library preparation protocol used for a unit. `STAR` can produce counts for unstranded (`none` - this is the default, e.g., Smart-seq2), forward oriented (`yes` e.g., QuantSeq) and reverse oriented (`reverse`) protocols. 
  - (optional, **but highly recommended**) metadata: additional sample metadata/annotation columns can/should be added and will be included on a sample-basis in the output sample annotation file (i.e., sequencing units from the same sample should have the same metadata and only differ in their `bam_file` column). Same rules as for `sample_name` apply.

Set workflow-specific `resources` or command line arguments (CLI) in the workflow profile `workflow/profiles/default.config.yaml`, which supersedes global Snakemake profiles.
