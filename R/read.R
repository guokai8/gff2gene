#' import gff3 or gtf file
#' @title read GFF file into Granges object
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors isSingleString
#' @param file filename for gtf or gff3
#' @param format character indicate the file format(gtf/gff3)
#' @return GRanges object
#' @author Kai Guo
#' @export
#'
importGFF <- function(file,format=c("auto","gtf","gff3")){
    format <- match.arg(format)
    if (format == "auto") {
        format <- .detect_file(file)
        if (!(format %in% c("gff3", "gff", "gtf")))
            stop(wmsg("Cannot detect whether 'file' is a GFF3 or GTF file. ",
                      "Please use the 'format' argument to specify the ",
                      "format (\"gff3\" or \"gtf\")."))
    }
    if (format == "gff3") {
        colnames <- GFF3_COLNAMES
    }
    else if (format == "gtf") {
        colnames <- GFF2_COLNAMES
    }
    else {
        colnames <- union(GFF3_COLNAMES, GFF2_COLNAMES)
    }
    message("Import genomic features from the file as a GRanges object ...\n ",
            appendLF = FALSE)
    gr <- import(file, format = format, colnames = colnames,
                 feature.type = GFF_FEATURE_TYPES)
    gr
}

#'
#' @title filter GRanges and add flanking regions
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @importFrom GenomicRanges GRanges
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors elementMetadata
#' @importFrom rlang `!!`
#' @importFrom rlang sym
#' @param gr GRanges object
#' @param gene gene name (set to NULL if specify the chromosome name and region)
#' @param chr chromosome name
#' @param region a vector include the start and the end location
#' @param exon a logical value indicating whether use exon information or not
#' @param cds a logical value indicating whether use cds information or not
#' @param gene_name a logical value indicating whether use
#'                  gene_name or not(for gtf)
#' @param format If not missing, should be one of “gff”,“gff3” or “gtf”.
#' @return GRangesInfo object include only select genes information
#' @export
#' @author Kai Guo
preGranges <- function(gr, gene = NULL,chr=NULL,region=NULL, exon = TRUE, cds = FALSE,
                       gene_name = FALSE,format = NULL)
{
    #gr <- .importGFF(file, type = "auto")
    if(is.null(format)){
        if("Parent" %in% names(mcols(gr))){
            format="gff3"
        }else{
            format="gtf"
        }
    }
    if(format == "gff3"){
        gene.col = "Parent"
    }else{
        if(isTRUE(gene_name)){
            gene.col = "gene_name"
        }else{
            gene.col = "gene_id"
        }
    }
    if(isTRUE(cds)){
        cols <- c("mRNA","transcript","CDS","three_prime_UTR","five_prime_UTR",
                  "start_codon","stop_codon")
    }else{
        cols <- c("mRNA","transcript","exon","three_prime_UTR","five_prime_UTR",
                  "start_codon","stop_codon")
    }
    cols <- intersect(cols,as.vector(gr$type))
    if(!is.null(chr)&!is.null(region)){
        gx <- GRanges(seqnames = chr, strand = "*", ranges = IRanges(start=region[1],end=region[2]))
        chr_r <- subsetByOverlaps(gff, gx, ignore.strand=TRUE)
        gene <- unique(sub('\\.\\d+','',unique(unlist(elementMetadata(chr_r)[,gene.col]))))
    }
    gr <- sapply(gene, function(x)gr %>% as.data.frame() %>%
                     filter(grepl(!!x,!!sym(gene.col))) %>%
                     filter(type %in% cols)%>%GRanges())%>%GRangesList()
    res <- new("GRangesInfo",
               gr = gr,
               gene = gene,
               format = format,
               mcol = gene.col,
               cols = cols)
    res
}
