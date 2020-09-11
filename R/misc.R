GFF3_COLNAMES <- c("source","type","phase","ID","Name","Parent","GeneID")
GFF2_COLNAMES<-c("source","type","phase",
                 "gene_id","gene_name","transcript_id","transcript_name")

#' copy from GenomicFeatures
.GENE_TYPES <- c("gene", "pseudogene", "transposable_element_gene")

.TX_TYPES <- c("transcript", "pseudogenic_transcript", "primary_transcript",
               "mRNA", "ncRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "tmRNA",
               "miRNA", "miRNA_primary_transcript",
               "RNase_P_RNA", "RNase_MRP_RNA", "SRP_RNA", "misc_RNA",
               "antisense_RNA", "antisense",
               "lnc_RNA", "antisense_lncRNA", "transcript_region",
               "pseudogenic_tRNA", "scRNA", "guide_RNA", "telomerase_RNA",
               "vault_RNA", "Y_RNA")

.EXON_TYPES <- c("exon", "pseudogenic_exon", "coding_exon",
                 "five_prime_coding_exon", "interior_coding_exon",
                 "three_prime_coding_exon", "exon_of_single_exon_gene",
                 "interior_exon", "noncoding_exon",
                 "five_prime_noncoding_exon", "three_prime_noncoding_exon")
.CDS_TYPES <- c("CDS", "transposable_element_CDS",
                "CDS_predicted", "edited_CDS")
.UTR_TYPES <- c("five_prime_UTR", "three_prime_UTR","UTR")

.STOP_CODON_TYPES <- c("start_codon","stop_codon")

### cory from GenomicFeatures
GFF_FEATURE_TYPES <- c(.GENE_TYPES, .TX_TYPES, .EXON_TYPES,
                       .CDS_TYPES, .STOP_CODON_TYPES,.UTR_TYPES)
GENE_FEATURE <- c(.EXON_TYPES,.CDS_TYPES,.UTR_TYPES,.STOP_CODON_TYPES)
#' detect file format for import gff or gff3 file
#' modified based on GenomicFeatures .detect_file_format function
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom rtracklayer FileForFormat
.detect_file <- function(file)
{
    if (isSingleString(file)) {
        file2 <- try(FileForFormat(file), silent=TRUE)
        if (inherits(file2, "try-error"))
            return(file_ext(file))
        file <- file2
    }
    if (is(file, "RTLFile")) {
        if (is(file, "GFF3File"))
            return("gff3")
        if (is(file, "GTFFile"))
            return("gtf")
        desc <- rtracklayer:::resourceDescription(file)
        if (is(file, "CompressedFile")) {
            ## Mimic what import,CompressedFile,missing,missing method does.
            desc <- file_path_sans_ext(desc)
        }
        type <- file_ext(desc)
        return(type)
    }
    stop(wmsg("Invalid 'file'. Must be a path to a file, or an URL, ",
              "or a GFF3File or GTFFile object."))
}
#' expand range from start and end
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom BiocGenerics which
#' @importFrom GenomicRanges `start<-`
#' @importFrom GenomicRanges `end<-`
#' @param gr GRanges object
#' @param upstream flanking length for upstream
#' @param downstream flanking length for downstream
#' @export
#' @author Kai Guo
expandRange = function(gr, upstream=1000, downstream=1000) {
    min = start(range(gr))
    max = end(range(gr))
    strand_is_minus = strand(gr) == "-"
    on_plus = which(!strand_is_minus)
    on_minus = which(strand_is_minus)
    start(gr)[on_plus & start(gr) == min] = start(gr)[on_plus &
                                                          start(gr) == min] - upstream
    end(gr)[on_minus & end(gr) == max] = end(gr)[on_minus &
                                                     end(gr) == max] + upstream
    end(gr)[on_plus & end(gr) == max] = end(gr)[on_plus & end(gr) == max] +
        downstream
    start(gr)[on_minus & start(gr) == min] = start(gr)[on_minus &
                                                           start(gr) == min] - downstream
    gr
}

#' unlist and add names
#' @param list a list include character vector
#' @author Kai Guo
.unlist.name <- function(list){
    list.name <- rep(names(list),times=unlist(lapply(list,length)))
    list <- unlist(list)
    names(list) <- list.name
    list
}
#' change mcols data colnames
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors `elementMetadata<-`
#' @importFrom S4Vectors elementMetadata
#' @param gr GRanges object
#' @param names new column names
#' @author Kai Guo
.colnames.mcol <- function(gr,names){
    #### need to be change
    names(elementMetadata(gr))[3:ncol(elementMetadata(gr))] <- names
    gr
}
