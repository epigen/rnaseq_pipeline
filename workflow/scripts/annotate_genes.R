library(tidyverse)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
# useful error messages upon aborting
library("cli")

#### config

# input
counts_path <- file.path(snakemake@input[["counts"]])
gtf_path <- file.path(snakemake@input[["gtf"]])
fasta_path <- file.path(snakemake@input[["fasta"]])

# output
gene_annot_path <- file.path(snakemake@output[["gene_annotation"]])

# params
species <- snakemake@params[["species"]]
version <- snakemake@params[["version"]]

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because 
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  mart <- tryCatch(
    {
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(species, "_gene_ensembl"),
          version = version,
        mirror = mart
      )
    },
    error = function(e) {
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        cli_abort(
          str_c(
            "Have tried all 4 available Ensembl biomaRt mirrors ",
            rounds,
            " times. You might have a connection problem, or no mirror is responsive.\n",
            "The last error message was:\n",
            message(e)
          )
        )
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}

# get quantified Ensembl gene IDs
counts <- read.table(counts_path, sep=',', header=1)

# annotate Ensembl gene IDs using biomaRt
gene_annot <- biomaRt::getBM(
            attributes = c( "ensembl_gene_id",
                            "version",
                            "source",
                            "external_gene_name",
                            "external_gene_source",
                            "description",
                            "gene_biotype"),
            filters = "ensembl_gene_id",
            values = counts$gene,
            mart = mart,
            )

#### get gc-content and gene length
# Code adapted from here: https://www.biostars.org/p/91218/#9612509
#Load the annotation and reduce it
GTF <- import.gff(gtf_path, format="gtf", feature.type="exon") # if we would want to inlcude introns use feature.type = c("exon","intron")
grl <- reduce(split(GTF, elementMetadata(GTF)$gene_id))
reducedGTF <- unlist(grl, use.names=T)
elementMetadata(reducedGTF)$gene_id <- rep(names(grl), elementNROWS(grl))

#Open the fasta file
FASTA <- FaFile(fasta_path)
open(FASTA)

#Add the GC numbers
elementMetadata(reducedGTF)$nGCs <- letterFrequency(getSeq(FASTA, reducedGTF), "GC")[,1]
elementMetadata(reducedGTF)$widths <- width(reducedGTF)

#Create a list of the ensembl_id/GC/length
calc_GC_length <- function(x) {
    nGCs = sum(elementMetadata(x)$nGCs)
    width = sum(elementMetadata(x)$widths)
    c(width, nGCs/width)
}
gc_length <- t(sapply(split(reducedGTF, elementMetadata(reducedGTF)$gene_id), calc_GC_length))
colnames(gc_length) <- c("exon_length", "exon_gc")
gc_length <- as.data.frame(gc_length)

#### merge and save gene annotations
gene_annot <- cbind(gene_annot, gc_length[gene_annot$ensembl_gene_id, ])
write.table(gene_annot, file=gene_annot_path, sep=",", quote=TRUE, row.names=FALSE)
