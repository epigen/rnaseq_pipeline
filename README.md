[![MrBiomics](https://img.shields.io/badge/MrBiomics-red)](https://github.com/epigen/MrBiomics/)
[![DOI](https://zenodo.org/badge/XXXXXX.svg)](https://zenodo.org/badge/latestdoi/XXXXXX)
[![](https://tokei.rs/b1/github/epigen/rnaseq_pipeline?category=code)]()
[![](https://tokei.rs/b1/github/epigen/rnaseq_pipeline?category=files)]()
[![GitHub license](https://img.shields.io/github/license/epigen/rnaseq_pipeline)](https://github.com/epigen/rnaseq_pipeline/blob/main/LICENSE)
![GitHub Release](https://img.shields.io/github/v/release/epigen/rnaseq_pipeline)
[![Snakemake](https://img.shields.io/badge/Snakemake->=8.25.3-green)](https://snakemake.readthedocs.io/en/stable/)

# RNA-seq Data Processing, Quantification & Annotation Pipeline
A [Snakemake 8](https://snakemake.readthedocs.io/en/stable/) workflow for end-to-end processing, quantification, and annotation of gene expression for RNA-seq experiments, starting from single- or paired-end reads within raw/unaligned/unmapped [uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format) files, including a comprehensive MultiQC report.

> [!NOTE]  
> This workflow adheres to the module specifications of [MrBiomics](https://github.com/epigen/MrBiomics), an effort to augment research by modularizing (biomedical) data science. For more details, instructions, and modules check out the project's repository.
>
> ⭐️ **Star and share modules you find valuable** 📤 - help others discover them, and guide our future work!

> [!IMPORTANT]  
> **If you use this workflow in a publication, please don't forget to give credit to the authors by citing it using this DOI [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) and acknowledging the [rna-seq-star-deseq2 workflow](https://github.com/snakemake-workflows/rna-seq-star-deseq2) DOI [10.5281/zenodo.4737358](https://doi.org/10.5281/zenodo.4737358) from which structure and code were adapted.**

![Workflow Rulegraph](./workflow/dags/rulegraph.svg)

# 🖋️ Authors
- [Stephan Reichl](https://github.com/sreichl)
- [Christoph Bock](https://github.com/chrbock)

# 💿 Software
This project wouldn't be possible without the following software and their dependencies.

| Software             | Reference (DOI / URL)                             |
| :-------------------: | :-----------------------------------------------: |
| Snakemake            | https://doi.org/10.12688/f1000research.29032.2    |
| STAR                 | https://doi.org/10.1093/bioinformatics/bts635|
| Samtools             | https://doi.org/10.1093/bioinformatics/btp352 |
| MultiQC              | https://doi.org/10.1093/bioinformatics/btw354     |
| RSeQC                | https://doi.org/10.1093/bioinformatics/bts356 |
| biomaRt              | https://doi.org/10.1038/nprot.2009.97 |
| EDASeq               | https://doi.org/10.1186/1471-2105-12-480 |
| gffutils             | https://github.com/daler/gffutils                 |
| pandas               | https://doi.org/10.5281/zenodo.3509134            |
| fastp                | https://doi.org/10.1093/bioinformatics/bty560   |
| pigz                 | https://zlib.net/pigz/                          |

# 🔬 Methods
This is a template for the Methods section of a scientific publication and is intended to serve as a starting point. Only retain paragraphs relevant to your analysis. References `[ref]` to the respective publications are curated in the software table above. Versions `(ver)` have to be read out from the respective conda environment specifications (`workflow/envs/*.yaml` file) or post-execution in the result directory (`rnaseq_pipeline/envs/*.yaml`). Parameters that have to be adapted depending on the data or workflow configurations are denoted in squared brackets e.g., `[X]`.

**Processing.** All unmapped BAM files for each sample were merged, converted to FASTQ format using Samtools [ref] (ver), processed for adapter trimming and quality filtering with fastp [ref] (ver) using parameters [`config["fastp_args"]`] (and adapters [`config["adapter_fasta"]`]), and finally de-interleaved (if paired-end) and compressed into separate R1/R2 files using shell commands and pigz [ref] (ver).

**Quantification.** Gene expression quantification was performed on the filtered and trimmed reads. The STAR aligner (ver) [ref] was utilized in `--quantMode GeneCounts` mode to count reads overlapping annotated genes based on the Ensembl [ver `config: ref: release`] gene annotation for the [`config: ref: species`] genome (build [`config: ref: build`]) and using parameters [`config["star_args"]`]. Library strandedness (`[none/yes/reverse]`, specified per sample in `config/annotation.csv`) was accounted for during counting. Counts for individual samples were aggregated into a single gene-by-sample count matrix using a custom Python script utilizing the pandas library (ver) [ref]. Quality control metrics from various tools, including fastp (ver) [ref], RSeQC (ver) [ref] and STAR (ver) [ref], were aggregated using MultiQC (ver) [ref].

**Annotation.** Gene annotations from Ensembl were retrieved using the R package biomaRt (ver) [ref]. Annotations included Ensembl gene ID, gene symbol (`external_gene_name`), gene biotype, and description. Additionally, exon-based GC content (`exon_gc`) and cumulative exon length (`exon_length`) were calculated for each gene using a custom R function adapted from EDASeq (ver) [ref], leveraging biomaRt (ver) [ref] to fetch exon coordinates and sequences. This exon-based approach was chosen as sequencing reads in poly(A)-selected libraries primarily derive from exonic regions, making these metrics more appropriate for downstream bias correction (e.g., Conditional Quantile Normalization - CQN) than whole-gene metrics. A sample annotation file was generated, integrating input annotations with QC metrics.

The processing and quantification described here was performed using a publicly available Snakemake [ver] (ref) workflow [[10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX)].

# 🚀 Features
This workflow offers several key advantages for RNA-seq analysis over existing pipelines:

*   **MrBiomics Module:** Designed for modularity and seamless integration with other [MrBiomics](https://github.com/epigen/MrBiomics/) analysis workflows (e.g., filtering/normalization, differential expression, unsupervised analysis) and includes example analysis recipes.
*   **Exon-Centric Annotation:** Calculates *exon-based* GC content and length, providing more accurate metrics for downstream bias correction (like CQN) for typical poly(A)-selected RNA-seq libraries.
*   **Robust Input Handling:** Starts directly from raw uBAM files and includes automated checks to ensure consistency between annotated read types (single/paired) and BAM file content.
*   **Efficiency Focused:** Optimized for performance and disk space using data streaming between processing steps and temporary intermediate files.
*   **User-Friendly:** Offers comprehensive documentation, clear configuration, standard Snakemake usage, and practical usage tips including detailed QC guidelines.
*   **Downstream-Ready Output:** Generates additionally to a count matrix, both gene and sample annotation files specifically formatted for easy use with common downstream analysis tools and modules.

---

The workflow performs the following steps that produce the outlined results:

- **Processing:**
    - Automatically verifies (using samtools view) that the `read_typ` (single/paired-end) specified in the annotation matches the actual flags within the input BAM files, preventing downstream errors (`.check_read_type/{sample}.done`).
    - Combines multiple input raw/unaligned/unmapped [uBAM](https://gatk.broadinstitute.org/hc/en-us/articles/360035532132-uBAM-Unmapped-BAM-Format) files per sample into a single stream (using `samtools merge`).
    - Converts the merged BAM stream into FASTQ format, handling paired-end interleaving (using `samtools fastq`).
    - Processes the FASTQ stream for adapter trimming and quality filtering using `fastp`, generating QC reports (`fastp/{sample}/`).
    > [!NOTE]  
    > `fastp` adapter auto-detection is disabled because we use STDIN mode (i.e., stream the data through pipes) to be disk space efficient.
    - De-interleaves the filtered FASTQ stream into separate compressed R1 and R2 files for paired-end data, or compresses directly for single-end data using shell commands and `pigz`.
- **Quantification:**
    - Uses STAR `GeneCounts` to quantify reads per gene based on the specified Ensembl reference genome and annotation (`star/{sample}/`).
    - Handles unstranded, forward-stranded, and reverse-stranded library protocols based on the `strandedness` column.
    - Aggregates counts into a single matrix (`counts/counts.csv`).
- **Annotation:**
    - Outputs gene annotations (`counts/gene_annotation.csv`).
      - Retrieves gene annotations (Ensembl ID, gene symbol, biotype, description) from Ensembl using `biomaRt`.
      - Calculates **exon-based** GC content and cumulative exon length for each gene, suitable for poly(A) selected libraries (e.g., .
      > [!NOTE]
      > Gene annotation can take a while since it depends on the availability of external data sources accessed via `biomaRt`.
      > GC-content and length are **exon-based**: In poly(A)‑selected libraries (such as Illumina TruSeq, Smart-seq or QuantSeq), the sequencing reads mainly come from exonic regions. Therefore, potential correction for GC bias and gene length should ideally use exon‑level GC content and effective exon length rather than whole‑gene metrics that include introns.
    - Outputs a sample annotation table containing sample-wise general MultiQC statistics (`counts/sample_annotation.csv`).
- **QC & Reporting:**
    - Employs RSeQC tools to generate key quality metrics like strand specificity and read distribution across genomic features (`rseqc/`).
    - Aggregates QC metrics from fastp, STAR and RSeQC into a single report using MultiQC (`report/multiqc_report.html`) with [AI summaries](https://seqera.io/blog/ai-summaries-multiqc/).

---

The workflow produces the following directory structure:
```
{result_path}/
└── rnaseq_pipeline/
    ├── report/
    │   └── multiqc_report.html         # Aggregated QC report for all samples
    ├── fastp/                          # fastp QC/filtering and adapter trimming outputs per sample
    ├── rseqc/                          # RSeQC output per sample
    ├── star/                           # STAR output per sample
    ├── counts/                         # Final quantification and annotation outputs
    │   ├── counts.csv                  # Aggregated gene count matrix (genes x samples)
    │   ├── sample_annotation.csv       # Annotation for samples (columns of counts.csv)
    │   └── gene_annotation.csv         # Annotation for genes (rows of counts.csv)
    ├── envs/                           # Exported Conda environment specifications
    └── configs/                        # Exported configuration and annotation files used for the run
```
> [!IMPORTANT]  
> Resources are downloaded automatically to `resources/{config::project_name}/rnaseq_pipeline/)`, are large (>`3GB`) and have to be manually removed if not needed anymore.


# 🛠️ Usage

Here are some tips for the usage of this workflow:

- Configure and run the workflow first for a few samples.
- Once everything works, run it for all samples.
- If you are unsure about the memory requirements for alignment (the most memory intensive step) provide the Snakemake parameter `--retries X`, where `X` denotes the number of retries, and with every retry the memory is increased to `attempts * config:mem`.
- To save disk space the intermediate gzipped FASTQ files are marked as temporary using Snakemake's `temp()` directive. To remove them upon successful completion you have to include the `--delete-temp-output` flag in your Snakemake command.

This workflow is written with Snakemake and its usage is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=epigen/rnaseq_pipeline).

# ⚙️ Configuration
Detailed specifications can be found here [./config/README.md](./config/README.md)

# 📖 Examples
Explore detailed examples showcasing module usage in comprehensive end-to-end analyses (including data, configuration, annotation and results) in our [MrBiomics Recipes](https://github.com/epigen/MrBiomics?tab=readme-ov-file#-recipes):
- [RNA-seq Analysis Recipe](https://github.com/epigen/MrBiomics/wiki/RNAseq-Analysis-Recipe)
  
# 🔍 Quality Control
Below are some guidelines for the manual quality control of each sample using the generated `MultiQC` report, but keep in mind that every experiment/dataset is different. Thresholds are general suggestions and may vary based on experiment type, organism, and library prep.

- **Alignment Rate (STAR):** % (Uniquely) Mapped Reads > 70-80%. Low rates might indicate contamination or reference issues.
- **Alignment Scores & Gene Counts (STAR):** High proportion of uniquely mapped reads assigned to exonic regions (e.g., >60-70% for poly(A) mRNA-seq). Low rates could suggest gDNA contamination or high intronic reads. In case of many intronic reads or non-poly(A) mRNA-seq protocols do not use exon-based gene annotations (gc-content and length).
- **Read Quality (fastp):** High average quality scores across reads after trimming. Ensure effective adapter removal.
- **Library Complexity (fastp/RSeQC):** % Duplication Rate should not be excessively high (highly variable, interpret in context of expression). Very high rates might indicate low input or PCR issues.
- **Strand Specificity (RSeQC):** For stranded protocols, >90-95% reads should match the expected strand.
- Inspect [**Genome Browser Tracks**](https://github.com/epigen/genome_tracks/) using UCSC Genome Browser (online) or IGV (local)
    - Compare all samples to the best, based on above's QC metrics.
    - Check cell type / experiment-specific markers or sex chromosome (`X`/`Y`) for expression as **positive controls**.
    - Check e.g., developmental regions for expression as **negative controls**.
- [**Unsupervised Analysis**](https://github.com/epigen/unsupervised_analysis) (e.g., PCA or UMAP)
    - Identify outliers/drivers of variation, especially in the control samples and within replicates.

# 🔗 Links
- [GitHub Repository](https://github.com/epigen/rnaseq_pipeline/)
- [GitHub Page](https://epigen.github.io/rnaseq_pipeline/)
- [Zenodo Repository](https://doi.org/10.5281/zenodo.XXXXXX)
- [Snakemake Workflow Catalog Entry](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/epigen%20rnaseq_pipeline.html)

# 📚 Resources
- Recommended compatible [MrBiomics Modules](https://github.com/epigen/MrBiomics/#-modules)
  - for upstream sample acquisition:
    - [Fetch Public Sequencing Data and Metadata Using iSeq](https://github.com/epigen/fetch_ngs/) to retrieve and prepare public data for downstream processing.
  - for downstream analyses:
      - [Genome Browser Track Visualization](https://github.com/epigen/genome_tracks/) for quality control and visual inspection/analysis of genomic regions/genes of interest or top hits.
      - [<ins>Sp</ins>lit, F<ins>ilter</ins>, Norma<ins>lize</ins> and <ins>Integrate</ins> Sequencing Data](https://github.com/epigen/spilterlize_integrate/) using the generated `counts.csv`, `sample_annotation.csv`, and `gene_annotation.csv` (especially the exon-based GC/length for CQN).
      - [Unsupervised Analysis](https://github.com/epigen/unsupervised_analysis) to explore sample relationships based on gene expression.
      - [Differential Analysis with limma](https://github.com/epigen/dea_limma) to identify differentially expressed genes between sample groups.
      - [Enrichment Analysis](https://github.com/epigen/enrichment_analysis) for biomedical interpretation of gene lists (e.g., from differential analysis).
- [fastp manual](https://github.com/OpenGene/fastp/blob/master/README.md)
- [STAR manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)

# 📑 Publications
The following publications successfully used this module for their analyses.
- [FirstAuthors et al. (202X) Journal Name - Paper Title.](https://doi.org/10.XXX/XXXX)
- ...

# ⭐ Star History

[![Star History Chart](https://api.star-history.com/svg?repos=epigen/rnaseq_pipeline&type=Date)](https://star-history.com/#epigen/rnaseq_pipeline&Date)

