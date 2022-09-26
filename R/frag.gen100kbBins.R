#' @title Generate 100kb bin fragment coverage dataframe
#'
#' @description
#' The 100kb bin fragment coverage dataframe including the following columns:
#'     "seqnames", "start", "end", "width", "strand", "Seqlength",
#'     "arm", "C.G", "short", "long", "ratio", "nfrags", "coverage",
#'     "short.corrected", "long.corrected", "nfrags.corrected", "ratio.corrected",
#'     "coverage.corrected", "combined", "short.predicted", "long.predicted",
#'     "nfrags.predicted", "ratio.predicted", "coverage.predicted",
#'     "mode", "mean", "median", "quantile.25", "quantile.75", "frag.gc",
#'     "id", and number of fragments with width of 90 up to 220
#'
#' @param object     fragSet object
#' @param filters    Character of path to filters rda file
#'                       (ENCFF356LFX.bed blacklist regions (910 ranges) to filter out)
#' @param gaps       Character of path to gaps rda file
#'                       (centromere, telomere and contig gaps (651 ranges) to filter out)
#' @param VNTRs      Character of path to VNTRs rda file
#'                       (variable number tandem repeats to (382 ranges) to filter out)
#' @param tiles      Character of path to tiles bed file
#'                       (30894 100kb bins template used to put fragments in for coverage calculation)
#' @param refGenome  Character of reference genome used
#'                       (default is "BSgenome.Hsapiens.UCSC.hg38")
#' @param chr.select Vector for selecting chromosomes to include in analysis
#'                       (default is paste0("chr",1:22))
#'
#' @return Returns a fragSet object with 100kb bin fragment coverage dataframe stored
#'     in "bins" slot, access with get100kbBins() function
#'
#' @example fset = gen100kbBins(fset,
#'                              filters=filters_path, gaps=gaps_path,
#'                              VNTRs=VNTRs_path, tiles=tiles_path,
#'                              refGenome="BSgenome.Hsapiens.UCSC.hg38",
#'                              chr.select=paste0("chr",1:22))
#'
#' @import tidyverse
#' @import GenomicRanges
#' @import IRanges
#' @importFrom utils read.table
#' @importFrom tibble as_tibble
#' @importFrom BSgenome getBSgenome
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom biovizBase GCcontent
#' @importFrom stats loess predict na.omit
#' @rdname gen100kbBins-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('gen100kbBins',
           function(object, filters, gaps, VNTRs, tiles,
                    refGenome="BSgenome.Hsapiens.UCSC.hg38",
                    chr.select=paste0("chr",1:22)) standardGeneric('gen100kbBins'))
#' @rdname gen100kbBins-methods
setMethod('gen100kbBins', 'fragSet',
          function(object,
                   filters, gaps, VNTRs, tiles,
                   refGenome="BSgenome.Hsapiens.UCSC.hg38",
                   chr.select=paste0("chr",1:22)){

            message("Loading filters, gaps, VNTRs and tiles (100kb bin template)")
            load(filters)
            load(gaps)
            load(VNTRs)
            hsapiens = getBSgenome(refGenome)
            AB = utils::read.table(tiles, col.names = c("chrom", "chromStart", "chromEnd", "Seqlength"))

            ## get 100kb bin template
            AB = makeGRangesFromDataFrame(AB, keep.extra.columns=TRUE)

            chrNums = as.numeric(gsub("chr","", chr.select))
            chromosomes = GRanges(chr.select,
                                  IRanges(0, seqlengths(hsapiens)[chrNums]))

            tcmeres = gaps.hg38[grepl("centromere|telomere", gaps.hg38$type)]

            ## removes warning: "GRanges object contains 24 out-of-bound ranges..."
            defaultW <- getOption("warn")
            options(warn = -1)
            arms = GenomicRanges::setdiff(chromosomes, tcmeres)
            options(warn = defaultW)
            arms = arms[-c(25,27,29,41,43)]

            armlevels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                          "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                          "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                          "19p", "19q","20p","20q","21q","22q")

            arms$arm = armlevels
            AB = AB[-queryHits(findOverlaps(AB, gaps.hg38))]
            AB = AB[queryHits(findOverlaps(AB, arms))]
            AB$arm = armlevels[subjectHits(findOverlaps(AB, arms))]

            seqinfo(AB) = seqinfo(hsapiens)[seqlevels(seqinfo(AB))]
            AB = trim(AB)
            AB$gc = GCcontent(hsapiens, AB)

            rm(chromosomes, tcmeres, arms, armlevels)

            message("Filter fragments into 100kb bins.")
            frags_Granges = getFragsGRanges(object)
            fragments = frags_Granges[-queryHits(findOverlaps(frags_Granges, filters.hg38))]
            fragments = fragments[queryHits(findOverlaps(fragments, AB))]
            fragments = fragments[-queryHits(findOverlaps(fragments, VNTRs.hg38))]
            rm(frags_Granges)

            message("Calculate coverage stats.")
            bamCov = coverage(fragments)
            mean = mean(bamCov)
            sd = sd(bamCov)
            max = max(bamCov)

            ## Set cutoff and report coverage
            cutoff = ceiling(median(quantile(bamCov, 0.9999)))
            adjust = ceiling(10*log2(median(mean)+1))
            cutoff = cutoff + ifelse(adjust > 5, adjust, 5)
            high_cov_regions = slice(bamCov, lower = cutoff)
            high_cov_regions = ranges(high_cov_regions)

            object@QC$high_cov_regions = high_cov_regions

            message("Generate coverage statistics.")
            ## (renamed from QC matrix)
            cov_stats = data.frame(mean, sd, max)

            object@QC$cov_stats = cov_stats
            rm(adjust, cutoff, max, mean, sd, high_cov_regions, bamCov, cov_stats)

            message("Count fragment sizes per bin...")
            w.all = width(fragments)
            fragments = fragments[which(w.all >= 90 & w.all <= 220)]
            w = width(fragments)
            frag.list = split(fragments, w)

            counts = sapply(frag.list, function(x) countOverlaps(AB, x))
            if(min(w) > 90) {
              m0 = matrix(0, ncol=min(w) - 90, nrow=nrow(counts),
                          dimnames=list(rownames(counts), 90:(min(w)-1)))
              counts = cbind(m0, counts)
            }

            if(max(w) < 220) {
              m1 = matrix(0, ncol=220 - max(w), nrow=nrow(counts),
                          dimnames=list(rownames(counts), (max(w)+1):220))
              counts = cbind(counts, m1)
            }

            olaps = findOverlaps(fragments, AB)
            bin.list = split(fragments[queryHits(olaps)], subjectHits(olaps))
            bingc = rep(NA, length(bin.list))
            bingc[unique(subjectHits(olaps))] = sapply(bin.list, function(x) mean(x$gc))
            gc = as.vector(AB$gc)
            while(length(bingc) < 26508){
              bingc = append(bingc, NA)
            }
            bingc = ifelse(is.na(bingc), gc, bingc)
            bingc = ifelse(bingc < min(gc), gc, bingc)
            bingc = ifelse(bingc > max(gc), gc, bingc)

            rm(filters.hg38, gaps.hg38, VNTRs.hg38, frag.list, fragments, bin.list, olaps)

            message("Calculate fragment statistics...")
            Mode = function(x) {
              ux = unique(x)
              ux[which.max(tabulate(match(x, ux)))]
            }
            modes = Mode(w)
            medians = median(w)
            q25 = quantile(w, 0.25)
            q75 = quantile(w, 0.75)

            short = rowSums(counts[,1:61])
            long = rowSums(counts[,62:121])
            ratio = short/long
            ratio[is.nan(ratio)] = NA
            ratio[is.infinite(ratio)] = NA
            nfrags = short+long
            coverage = nfrags/sum(nfrags, na.rm=TRUE)

            message("GC correction and prediction...")
            short.corrected = gc.correct(short, bingc)
            long.corrected = gc.correct(long, bingc)
            short.predicted = gc.pred(short, bingc)
            long.predicted = gc.pred(long, bingc)
            nfrags.predicted = gc.pred(short+long, bingc)
            coverage.predicted = gc.pred(coverage, bingc)

            short.corrected = ifelse(short.corrected <= 0, 0, short.corrected)
            long.corrected = ifelse(long.corrected <= 0, NA, long.corrected)
            nfrags.corrected = short.corrected+long.corrected
            ratio.corrected = short.corrected/long.corrected
            coverage.corrected = nfrags.corrected/sum(nfrags.corrected, na.rm=TRUE)
            combined = ratio.corrected*coverage.corrected

            ## Append fragment information
            AB$short = short
            AB$long = long
            AB$ratio = ratio
            AB$nfrags = short+long
            AB$coverage = coverage
            AB$short.corrected = short.corrected
            AB$long.corrected = long.corrected
            AB$nfrags.corrected = nfrags.corrected
            AB$ratio.corrected = ratio.corrected
            AB$coverage.corrected = coverage.corrected
            AB$combined = combined
            AB$short.predicted = short.predicted
            AB$long.predicted = long.predicted
            AB$nfrags.predicted = nfrags.predicted
            AB$ratio.predicted = short.predicted/long.predicted
            AB$coverage.predicted = coverage.predicted

            AB$mode = modes
            AB$mean = round(mean(w), 2)
            AB$median = medians
            AB$quantile.25 = q25
            AB$quantile.75 = q75
            AB$frag.gc = bingc

            sampleID=getSampleID(object)
            AB$id = sampleID

            for(i in 1:ncol(counts)){
              elementMetadata(AB)[,colnames(counts)[i]] = counts[,i]
            }

            rm(coverage, coverage.corrected, coverage.predicted, gc, i, long, long.corrected, long.predicted,
               medians, modes, nfrags, nfrags.corrected, nfrags.predicted, q25, q75, ratio, ratio.corrected,
               short, short.corrected, short.predicted, Mode, combined, bingc, w, counts)

            message("Generate raw 100kb bin object.")
            tib.list = as_tibble(AB)

            object@bins$kb100bins = tib.list
            rm(AB)

            return(object)
          })

## GC correct helper function
#' @rdname gen100kbBins-methods
#' @export
gc.correct = function(coverage, bias) {
  i = seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend = loess(coverage ~ bias, na.action = na.omit)
  coverage.model = loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred = predict(coverage.model, bias)
  coverage.corrected = coverage - coverage.pred + median(coverage, na.rm=TRUE)
}

## GC prediction helper function
#' @rdname gen100kbBins-methods
#' @export
gc.pred = function(coverage, bias) {
  i = seq(min(bias, na.rm=TRUE), max(bias, na.rm=TRUE), by = 0.001)
  coverage.trend = loess(coverage ~ bias, na.action = na.omit)
  coverage.model = loess(predict(coverage.trend, i) ~ i, na.action = na.omit)
  coverage.pred = predict(coverage.model, bias)
}

## loading extdata
# .onLoad = function (libname, pkgname) {
#   filters = system.file("extdata", "filters.hg38.rda", package = "fragmentomicsBEDPE")
#   assign('filters', filters, envir = topenv())
#
#   gaps = system.file("extdata", "gaps.hg38.rda", package = "fragmentomicsBEDPE")
#   assign('gaps', gaps, envir = topenv())
#
#   VNTRs = system.file("extdata", "VNTRs.hg38.rda", package = "fragmentomicsBEDPE")
#   assign('VNTRs', VNTRs, envir = topenv())
#
#   tiles = system.file("extdata", "hg38_tiles.bed", package = "fragmentomicsBEDPE")
#   assign('tiles', tiles, envir = topenv())
#
#   healthy = system.file("extdata", "healthy.median.hg38.rda", package = "fragmentomicsBEDPE")
#   assign('healthy', healthy, envir = topenv())
# }
#
# .onLoad()
