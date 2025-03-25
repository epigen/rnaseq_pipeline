library(biomaRt)
library(tidyverse)
library(EDASeq)
# useful error messages upon aborting
library("cli")

# adapted from EDASeq::getGeneLengthAndGCContent to work with specific ensembl
# otherwise data does not match or genes are not found
# returns exon gc content and exon lengths of transcripts for e.g., conditional quantile normalization
getGeneLengthAndGCContent <- function(id, id.type, ensembl){
    inp.id <- id
    message( paste0( "Downloading sequence",
        ifelse(length(id) > 1, "s", ""), " ..."))
    if(length(id) > 100) message("This may take a few minutes ...")

    # download sequence
    # (1) get exon coordinates
    attrs <- c(id.type, "ensembl_exon_id",
        "chromosome_name", "exon_chrom_start", "exon_chrom_end")
    coords <- getBM(filters=id.type, attributes=attrs, values=id, mart=ensembl)
    id <- unique(coords[,id.type])
    coords <- GRangesList(sapply(id,
        function(i)
        {
            i.coords <- coords[coords[,1]== i, 3:5]
            g <- GRanges(i.coords[,1], IRanges(i.coords[,2],i.coords[,3]))
            return(g)
        }), compress=FALSE)
    coords <- reduce(coords)
    len <- sum(width(coords))

    # (2) get genes and sequences
    sel <- c(id.type, "start_position", "end_position")
    gene.pos <- getBM(attributes = sel, filters=id.type, values=id,
                      mart=ensembl)
    gene.seqs <- getSequence(id=id,
        type=id.type, seqType="gene_exon_intron", mart=ensembl)

    # (3) get exonic sequences and correspondig GC content
    gc.cont <- sapply(id,
        function(i)
        {
            # exon coordinates, gene position & sequence for current id i
            ecoords <- coords[[i]]
            gpos <- gene.pos[gene.pos[,id.type] == i,
                    c("start_position", "end_position")]
            gseq <- DNAString(
                gene.seqs[gene.seqs[,id.type] == i, "gene_exon_intron"])

            # exon coordinates relative to gene position
            start <- start(ranges(ecoords)) - gpos[1,1] + 1
            end <- end(ranges(ecoords)) - gpos[1,1] + 1
            eseq <- gseq[IRanges(start, end)]
            gc.cont <- sum(alphabetFrequency(eseq, as.prob=TRUE)[c("C","G")])
            return(gc.cont)
        }
    )

    res <- cbind(len, gc.cont)
    colnames(res) <- c("length", "gc")
    rownames(res) <- id

    # (4) order according to input ids
    not.found <- !(inp.id %in% rownames(res))
    na.col <- rep(NA, sum(not.found))
    rn <- c(rownames(res), inp.id[not.found])
    res <- rbind(res, cbind(na.col, na.col))
    rownames(res) <- rn
    res <- res[inp.id,]
    return(res)
}

#### config

# input
counts_path <- file.path(snakemake@input[["counts"]])

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
df <- read.table(counts_path, sep=',', header=1)

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
            values = df$gene,
            mart = mart,
            )

# get gc-content and gene length (this takes some time)
gc_length <- getGeneLengthAndGCContent(
    id = gene_annot$ensembl_gene_id,
    id.type = 'ensembl_gene_id',
    ensembl = mart
)

# merge and save gene annotations
gene_annot <- cbind(gene_annot, gc_length[gene_annot$ensembl_gene_id, ])
write.table(gene_annot, file=gene_annot_path, sep=",", quote=TRUE, row.names=FALSE)
