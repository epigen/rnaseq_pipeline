# RNA-seq Data Processing, Quantification & Annotation Pipeline


Note: fastp adapter auto-detection is disabled because we use STDIN mode (i.e., stream the data through pipes) to be disk space efficient.

Note: Gene annotation takes a while as it depends on external data sources via accessed via biomaRt

Note: GC-content and length are exon based(!)
In poly(A)‑selected libraries, the sequencing reads mainly come from exonic regions.
Therefore, correcting for GC bias (and gene length) should ideally use exon‑level GC content and effective exon length rather than whole‑gene metrics that include introns.

Protocols focusing on poly(A) mRNA (poly‑A selection) according to o3-mini-high:

Illumina TruSeq Stranded mRNA Library Prep: Enriches for polyadenylated transcripts via oligo‑dT selection. 
BMCGENOMICS.BIOMEDCENTRAL.COM

NEBNext Poly(A) mRNA Isolation-based Kits: Typically use magnetic bead‑based poly(A) capture prior to library construction.

SMART‑Seq (e.g. SMART‑Seq2, SMART‑Seq v4): Full‑length cDNA synthesis from poly(A)‑tailed mRNAs, widely used for low‑input and single‑cell RNA‑seq.

Poly(A)‑ClickSeq: A variant of ClickSeq that selectively captures poly(A)‑transcripts without fragmentation.

10X Genomics Chromium 3′ Single‑Cell RNA‑Seq: Relies on poly(A) priming for cell barcoding and mRNA capture.

QuantSeq, developed by Lexogen, is a 3′ mRNA‑sequencing protocol that uses oligo‑dT priming to capture poly(A)‑tailed mRNAs.