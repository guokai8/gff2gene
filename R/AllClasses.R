##' Class "GRangesInfo"
##' This class represents the extract gene infomation from
##' @slot gr GrangesList of gff or gff3 files with filter genes
##' @slot gene vector of character of genes
##' @slot gene Gene IDs
##' @slot format a character gff/gtf or gff3
##' @slot mcol a character specify the gene id column
##' @slot cols a vector include gene structure information
##' @exportClass GRangesInfo
##' @author Kai Guo
##' @keywords classes
setClass("GRangesInfo",
         representation=representation(
             gr         = "GRangesList",
             gene       = "character",
             format     = "character",
             mcol       = "character",
             cols       = "character"
         )
)
