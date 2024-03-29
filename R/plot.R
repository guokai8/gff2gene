#' @title plot Gene Structure with gff
#' @importFrom IRanges subsetByOverlaps
#' @importFrom IRanges IRanges
#' @importFrom BiocGenerics start
#' @importFrom BiocGenerics end
#' @importFrom BiocGenerics strand
#' @importFrom BiocGenerics Reduce
#' @importFrom GenomicRanges mcols
#' @importFrom grid unit
#' @importFrom grid grid.rect
#' @importFrom grid grid.lines
#' @importFrom grid popViewport
#' @importFrom grid pushViewport
#' @importFrom grid grid.newpage
#' @importFrom grid grid.points
#' @importFrom grid grid.text
#' @importFrom grid grid.raster
#' @importFrom grid viewport
#' @importFrom grid convertHeight
#' @importFrom grid grid.legend
#' @importFrom grid grid.layout
#' @importFrom grid grid.yaxis
#' @importFrom grid grid.xaxis
#' @importFrom grid convertWidth
#' @importFrom grid convertHeight
#' @importFrom grid stringWidth
#' @importFrom grid arrow
#' @importFrom grid gpar
#' @importFrom grid grid.pretty
#' @importFrom grid grid.clip
#' @importFrom GenomicRanges GRangesList
#' @importFrom S4Vectors `elementMetadata<-`
#' @importFrom S4Vectors elementMetadata
#' @importFrom IRanges disjoin
#' @importFrom IRanges width
#' @param obj Granges object include gff file
#' @param dat Granges object include GWAS results
#' @param color color for the point
#' @param chr chromosome name
#' @param region a vector include the start and the end location
#' @param gene gene name
#' @param upstream upstream bp
#' @param downstream downstream bp
#' @param exon exon color
#' @param intron intron color
#' @param utr3 3'utr color
#' @param utr5 5'utr color
#' @param gene.label.size size of the gene label
#' @param arrow.col arrow color
#' @param arrow.fill arrow fill
#' @export
#' @author Kai Guo
geneplot <- function(obj, id = "gene_id", color=NULL, chr = NULL, region = NULL, gene = NULL,
                     geneOnly = FALSE, exonOnly = FALSE,
                     upstream = 1000, downstream = 1000, 
                     exon = "darkgreen", utr3 = "cyan4", utr5 ="cyan4",intron = "black", gene.label.size = 0.5,
                     ylab="Gene", legend.lab = NULL,
                     ylab.size = 0.9, arrow.col ="lightblue",arrow.fill = "lightblue"){
    #extrack object
    if(length(intersect(id,GFF3_COLNAMES))==0){
        id <- intersect(id,GFF2_COLNAMES)
    }else{
        id <- "Parent"
    }
    track <- obj[,c("type",id)]
    names(elementMetadata(track))<-c("feature","featureID")
    track <- track[track$feature%in%GENE_FEATURE,]
    track$featureID<-unlist(track$featureID)
    track$feature<-as.factor(as.character(track$feature))
   # names(elementMetadata(dat))[1]<-"scores"
    unit <- 0.25
    y <- 0.5
    ## whether only draw some genes
    if(!is.null(gene)){
        if(length(gene)>1){
            gene = paste(gene, collapse="|")
        }
        track <- track[grepl(gene,track$featureID,ignore.case = T),]
    }
    chr_r <- range(track,ignore.strand=T)
    ## draw only region
    if(!(is.null(chr) & is.null(region))){
        gr <- GRanges(seqnames = chr, strand = "*", ranges = IRanges(start=region[1],end=region[2]))
        chr_r <- subsetByOverlaps(gr, chr_r, ignore.strand=TRUE)
        track <- subsetByOverlaps(track, chr_r,ignore.strand=TRUE)
    }else{
        chr_r <- expandRange(chr_r, upstream = upstream,downstream = downstream)
    }
    ## start and end
    start <- start(chr_r)
    end <- end(chr_r)
    xscale <- c(start,end)
    # yscale
    ### layout
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(2,3,
                                                   widths =unit(c(0.5,8, 1.5),
                                                                units = c("null","null","null")),
                                                   heights = unit(c(9.5,0.5),units =c("null","null","null")))))
    #####################################################
    ## gene structure
    ## modify from trackViewer plotGeneTrack function
        vp <- viewport(layout.pos.row = 1,layout.pos.col = 2, xscale = xscale, y=0, height = 1,width = 0.8)
        pushViewport(vp)
        vp2 <- viewport(xscale = xscale, y=unit(0.6, "npc"),height = unit(0.8, "npc"))
        pushViewport(vp2)
        trs <- split(track, track$featureID)
        rgs <- unlist(GRangesList(sapply(trs, range)))
        if(length(trs)!=length(rgs)){
            stop("some genes are splited into multiple regions, eg multple strand.")
        }
        strand <- as.character(strand(rgs))
        dj <- disjoin(rgs, with.revmap = TRUE)
        lineId <- data.frame(id = seq_along(rgs), line = 0)
        revmap <- dj$revmap
        for(i in seq_along(revmap)){
            revm <- sort(revmap[[i]])
            if(lineId[revm[1], "line"] == 0){
                lineId[revm[1], "line"] <- 1
            }
            for(j in seq_along(revm)[-1]){
                lineId[revm[j], "line"] <- lineId[revm[j-1], "line"] + 1
            }
        }
        totalLines <- max(lineId$line)
        ## check height of plot region
        unity <- convertHeight(unit(1, "npc"), unitTo = "lines", valueOnly = TRUE)
        doLabels <- unity >= totalLines
        eachLineHeight <- 1/totalLines
        currLineBottom <- 1
        for(i in seq.int(totalLines)){
            vp <- viewport(y = currLineBottom - eachLineHeight/2,
                           height = eachLineHeight, xscale = xscale)
            pushViewport(vp)
            if(doLabels){
                gene_y <- .75
                gene_h <- .5
            }else{
                gene_y <- .5
                gene_h <- 1
            }
            trs.sub <- trs[which(lineId$line == i)]
            rgs.sub <- rgs[which(lineId$line == i)]
            str.sub <- strand[which(lineId$line == i)]
            oid <- order(start(rgs.sub))
            trs.sub <- trs.sub[oid]
            rgs.sub <- rgs.sub[oid]
            str.sub <- str.sub[oid]
            stringStopPos <- c(0, 0)
            for(j in seq_along(trs.sub)){
                curr_trs <- trs.sub[[j]]
                curr_rg <- rgs.sub[j]
                curr_str <- str.sub[j]
                str_neg <- curr_str == "-"
                ## plot centr line
                grid.lines(x = c(start(curr_rg), end(curr_rg)),
                           y = c(gene_y, gene_y),
                           gp = gpar(col = intron),
                           default.units = "native")
                exons <- curr_trs[curr_trs$feature%in%.EXON_TYPES]
                ## plot exons by the feature type
                grid.rect(x = (start(exons) + end(exons))/2, y = gene_y,
                          width = width(exons),
                          height = gene_h/ifelse(exons$feature %in%.EXON_TYPES , 2, 4),
                          gp = gpar(col = NA, fill = exon),
                          default.units = "native")
                utrs3 <- curr_trs[curr_trs$feature%in%c("utr3", "three_prime_UTR")]
                if(length(utrs3) > 0){
                    grid.rect(x = (start(utrs3) + end(utrs3))/2, y = gene_y,
                              width = width(utrs3),
                              height = gene_h/ifelse(utrs3$feature %in% c("utr3","three_prime_UTR"), 2, 4),
                              gp=gpar(col = NA, fill = utr3),
                              default.units = "native")
                }
                # utr5
                utrs5 <- curr_trs[curr_trs$feature%in%c("utr5", "five_prime_UTR")]
                if(length(utrs5) > 0){
                    grid.rect(x = (start(utrs5) + end(utrs5))/2, y = gene_y,
                              width = width(utrs5),
                              height = gene_h/ifelse(utrs5$feature %in% c("utr5", "five_prime_UTR"), 2, 4),
                              gp=gpar(col = NA, fill = utr5),
                              default.units = "native")
                }
                utrs <- curr_trs[curr_trs$feature%in%c("UTR")]
                if(length(utrs) > 0){
                    grid.rect(x = (start(utrs) + end(utrs))/2, y = gene_y,
                              width = width(utrs),
                              height = gene_h/ifelse(utrs$feature %in% c("UTR"), 2, 4),
                              gp=gpar(col = NA, fill = "cyan4"),
                              default.units = "native")
                }
                ## plot direction at TSS
                stringW <- convertWidth(stringWidth(names(curr_rg)), unitTo = "native", valueOnly = TRUE)
                pushViewport(viewport(x = ifelse(!str_neg, start(curr_rg), end(curr_rg)),
                                      y = gene_y,
                                      width = unit(max(min(gene_h/2, 2*width(curr_rg)/abs(diff(xscale))),
                                                       convertWidth(unit(8, "lines"),
                                                                    unitTo = "npc",
                                                                    valueOnly = TRUE)), "snpc"),
                                      height = unit(gene_h/2, "snpc"),
                                      just = c(ifelse(!str_neg, 0, 1), 1),
                                      default.units = "native"))
                if(!str_neg){
                    grid.lines(x = unit(c(0, 0, 1), "npc"),
                               y = unit(c(.5, 0, 0), "npc"),
                               arrow = arrow(type="closed", angle = 15, length = unit(.2, "lines")),
                               gp = gpar(col = arrow.col, fill = arrow.fill))
                    ## add gene name at TSS
                    if(doLabels){
                        if(start(curr_rg) >= stringStopPos[1]){
                            # grid.text(label = names(curr_rg), x = 0, y = 0, hjust = 0, vjust=1.5,
                            #              gp = gpar(track@style@ylabgp))
                            grid.text(label = names(curr_rg), x = 0, y = 0, hjust = 0, vjust = 1.5,
                                      gp = gpar(cex = gene.label.size,col = "black"), check.overlap = TRUE)
                            stringStopPos[1] <- start(curr_rg) + stringW
                        }else{
                            if(start(curr_rg) >= stringStopPos[2]){
                                grid.text(label = names(curr_rg), x = 0, y = 1, hjust = 0, vjust = -1,
                                          gp = gpar(cex = gene.label.size,col="black"), check.overlap = TRUE)
                                stringStopPos[2] <- start(curr_rg) + stringW
                            }
                        }
                    }
                }else{
                    grid.lines(x = unit(c(1, 1, 0), "npc"),
                               y = unit(c(.5, 0, 0), "npc"),
                               arrow = arrow(type = "closed", angle = 15, length = unit(.2, "lines")),
                               gp=gpar(col = arrow.col, fill = arrow.col))
                    ## add gene name at TSS
                    if(doLabels){
                        if(end(curr_rg) - stringW >= stringStopPos[1]){
                            grid.text(label = names(curr_rg), x = 1, y = 0, hjust = 1, vjust=1.5,
                                      gp = gpar(cex=gene.label.size, col="black"), check.overlap = TRUE)
                            stringStopPos[1] <- end(curr_rg)
                        }else{
                            if(end(curr_rg)-stringW >= stringStopPos[2]){
                                grid.text(label = names(curr_rg), x = 1, y = 1, hjust = 1, vjust=-1,
                                          gp = gpar(cex = gene.label.size, col = "black"), check.overlap = TRUE)
                                stringStopPos[2] <- end(curr_rg)
                            }
                        }
                    }
                }
                popViewport()
            }
            popViewport()
            currLineBottom <- currLineBottom - eachLineHeight
        }
        # xtick
        len <- end - start + 1
        if(len > 50000){
            xticklab <- paste0(round(grid.pretty(range = c(start, end))/1000000, 2), "Mb")
        }else if(len > 5000 & len <= 50000){
            xticklab <- paste0(round(grid.pretty(range = c(start, end))/1000, 2), "Kb")
        }else{
            xticklab <- grid.pretty(range = c(start, end))
        }
        grid.xaxis(at = grid.pretty(range = c(start,end)),label = xticklab)
        popViewport()
        popViewport()
        vp <- viewport(layout.pos.row = 1, layout.pos.col = 3, x=0.5, y=0.5, height = 0.5, width = 0.8)
        pushViewport(vp)
        vp3 <- viewport(x = unit(0.5, "npc"), y=unit(0.5, "npc"),height = unit(0.8, "npc"))
        pushViewport(vp3)
        if(length(utrs) > 0){
            grid.legend(labels=c("UTR","Intron","Exon"),ncol=1,byrow=TRUE, vgap=unit(.3, "lines"),
                        hgap=unit(.3, "lines"), pch=22, gp=gpar(col="white", fill=c("cyan4",intron,exon)))
        }else{
            grid.legend(labels=c("5'-UTR","3'-UTR","Intron","Exon"),ncol=1,byrow=TRUE, vgap=unit(.3, "lines"),
                        hgap=unit(.3, "lines"), pch=22, gp=gpar(col="white", fill=c(utr5,utr3,intron,exon)))
        }
        popViewport()
        popViewport()
    grid.clip()
    return(invisible())
}
