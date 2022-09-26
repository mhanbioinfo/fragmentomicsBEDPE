#' @title Write out dataframe results as text files
#'
#' @description
#' Write out mitochondrial report, high coverage regions, coverage statistics,
#'     100kb bin, 5Mb bin, 5Mb distance to healthy median, 5Mb summary files.
#'
#' @param object fragSet object
#' @param outdir Character of full path to output directory
#'
#' @return None
#'
#' @example outputAllResults(fset, outdir=outdir_path)
#'
#' @import tidyverse
#' @importFrom rtracklayer export.bed
#' @importFrom utils write.table
#' @rdname outputAllResults-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('outputAllResults',
           function(object, outdir) standardGeneric('outputAllResults'))
#' @rdname outputAllResults-methods
setMethod('outputAllResults', 'fragSet',
          function(object, outdir){
            dir.create(outdir, showWarnings = F)
            sampleID = getSampleID(object)
            message(paste0("Writing out all results for: ", sampleID, " ."))
            write.table(getMitoReport(object),
                        file.path(outdir, paste0(sampleID, "_result_mito_stats.txt")),
                        sep = "\t", quote = F, row.names = F)
            rtracklayer::export.bed(getHighCovRegions(object),
                                    file.path(outdir, paste0(sampleID, "_qc_highcovregions.bed")))
            write.table(getCovStats(object),
                        file.path(outdir, paste0(sampleID, "_qc_cov_stats.txt")),
                        sep = "\t", quote = F)
            write.table(get100kbBins(object),
                        file.path(outdir, paste0(sampleID, "_result_100kb_bins.txt")),
                        sep = "\t", quote = F, row.names = F)
            write.table(get5MbBins(object),
                        file.path(outdir, paste0(sampleID, "_result_5Mb_bins.txt")),
                        sep = "\t", quote = F, row.names = F)
            write.table(get5MbDist(object),
                        file.path(outdir, paste0(sampleID, "_result_5Mb_dist.txt")),
                        sep = "\t", quote = F, row.names = F)
            write.table(get5MbSummary(object),
                        file.path(outdir, paste0(sampleID, "_result_5Mb_summary.txt")),
                        sep = "\t", quote = F, row.names = F)
          })

#' @title Write out plots related to GC bias correction
#'
#' @description
#' Write out plots before and after GC bias correction
#'
#' @param object fragSet object
#' @param outdir Character of full path to output directory
#'
#' @return None
#'
#' @example outputGCcorrectPlots(fset, outdir=outdir_path)
#'
#' @import tidyverse
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics par smoothScatter
#' @rdname outputGCcorrectPlots-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('outputGCcorrectPlots',
           function(object, outdir) standardGeneric('outputGCcorrectPlots'))
#' @rdname outputGCcorrectPlots-methods
setMethod('outputGCcorrectPlots', 'fragSet',
          function(object, outdir){
            dir.create(outdir, showWarnings = F)
            sampleID = getSampleID(object)
            tib.list = dplyr::select(object@bins$kb100bins, -matches("X"))

            message("Plot GC Correction metrics.")
            pdf(file = file.path(outdir, paste0(sampleID, "_plot_GC_metrics.pdf")))
            par(mfrow=c(2,2))
            ## short
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$short,
                          main = "short",
                          xlab = "frag_GC",
                          ylab = "short")
            ## short corrected
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$short.corrected,
                          main = "short corrected",
                          xlab = "frag_GC",
                          ylab = "short_corrected")
            ## short vs short predicted
            smoothScatter(x = tib.list$short.predicted,
                          y = tib.list$short,
                          main = "short predicted vs actual",
                          xlab = "short_predicted",
                          ylab = "short")
            ## short corrected vs short predicted
            smoothScatter(x = tib.list$short.predicted,
                          y = tib.list$short.corrected,
                          main = "short predicted vs corrected",
                          xlab = "short_predicted",
                          ylab = "short_corrected")
            ## long
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$long,
                          main = "long",
                          xlab = "frag_GC",
                          ylab = "long")
            ## long corrected
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$long.corrected,
                          main = "corrected long",
                          xlab = "frag_GC",
                          ylab = "long_corrected")
            ## long vs long predicted
            smoothScatter(x = tib.list$long.predicted,
                          y = tib.list$long,
                          main = "long predicted vs actual",
                          xlab = "long_predicted",
                          ylab = "long")
            ## long corrected vs long predicted
            smoothScatter(x = tib.list$long.predicted,
                          y = tib.list$long.corrected,
                          main = "long predicted vs corrected",
                          xlab = "long_predicted",
                          ylab = "long_corrected")
            ## total fragments
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$nfrags,
                          main = "nfrags",
                          xlab = "frag_GC",
                          ylab = "nfrags")
            ## corrected total fragments
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$nfrags.corrected,
                          main = "corrected nfrags",
                          xlab = "frag_GC",
                          ylab = "nfrags_corrected")
            ## fragments vs predicted fragments
            smoothScatter(x = tib.list$nfrags.predicted,
                          y = tib.list$nfrags,
                          main = "nfrags predicted vs actual",
                          xlab = "nfrags_predicted",
                          ylab = "nfrags")
            ## corrected fragments vs predicted fragments
            smoothScatter(x = tib.list$nfrags.predicted,
                          y = tib.list$nfrags.corrected,
                          main = "nfrags predicted vs corrected",
                          xlab = "nfrags_predicted",
                          ylab = "nfrags_corrected")
            ## ratios
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$ratio,
                          main = "ratios",
                          xlab = "frag_gc",
                          ylab = "ratio")
            ## corrected ratios
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$ratio.corrected,
                          main = "corrected ratios",
                          xlab = "frag_gc",
                          ylab = "ratio_corrected")
            ## predicted ratios
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$ratio.predicted,
                          main = "predicted ratios",
                          xlab = "frag_gc",
                          ylab = "ratio_predicted")
            ## bin GC vs frag GC
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$C.G,
                          main = "GC content",
                          xlab = "frag_GC",
                          ylab = "bin_GC")
            ## coverage
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$coverage,
                          main = "coverage",
                          xlab = "frag_GC",
                          ylab = "coverage")
            ## corrected coverage
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$coverage.corrected,
                          main = "corrected coverage",
                          xlab = "frag_GC",
                          ylab = "coverage_corrected")
            ## predicted coverage
            smoothScatter(x = tib.list$frag.gc,
                          y = tib.list$coverage.predicted,
                          main = "predicted coverage",
                          xlab = "frag_GC",
                          ylab = "coverage_predicted")
            dev.off()
          })

#' @title Write out plots related to fragmentation profile
#'
#' @description
#' Write out plots including ...fragment.pdf (short/long fragment ratio in 5Mb bins),
#'     ...coverage.pdf (short fragment coverage in 5Mb bins),
#'     ...coverage.pdf (ratio*coverage in 5Mb bins).
#'
#' @param object fragSet object
#' @param healthy Character of path to healthy median rda file
#' @param outdir Character of full path to output directory
#'
#' @return None
#'
#' @example outputFragPlots(fset, healthy=healthy_path, outdir=outdir_path)
#'
#' @import tidyverse
#' @import ggplot2
#' @importFrom stats setNames
#' @rdname outputFragPlots-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('outputFragPlots',
           function(object, healthy, outdir) standardGeneric('outputFragPlots'))
#' @rdname outputFragPlots-methods
setMethod('outputFragPlots', 'fragSet',
          function(object, healthy, outdir){
            dir.create(outdir, showWarnings = F)
            load(healthy)
            id = getSampleID(object)
            df.fr3 = object@bins$mb5bins

            ## Set themes and plot layouts
            mytheme = theme_classic(base_size=12) + theme(
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              strip.text.x=element_text(size=11),
              strip.text.y=element_text(size=12),
              axis.title.x=element_text(face="bold", size=17),
              axis.title.y=element_text(size=15),
              axis.text.y=element_text(size=15),
              plot.title=element_text(size=15),
              legend.position="none",
              legend.title=element_text(size=10),
              legend.text=element_text(size=10),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              strip.background=element_rect(fill="white", color="white"),
              panel.spacing.x=unit(0.1, "lines"))

            armlevels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                          "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                          "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                          "19p", "19q","20p","20q","21q","22q")
            df.fr3$arm = factor(df.fr3$arm, levels=armlevels)
            healthy_median$arm = factor(healthy_median$arm, levels=armlevels)

            arm = df.fr3 %>% group_by(arm) %>%
              dplyr::summarize(n=n()) %>%
              mutate(arm = as.character(arm))
            small.arms = setNames(c("", "10q", "", "12q", "", "16",
                                    "", "17q", "", "18q",
                                    "", "", "", "",
                                    "", ""),
                                  c("10p", "10q", "12p", "12q", "16p", "16q",
                                    "17p", "17q", "18p", "18q",
                                    "19p", "19q", "20p", "20q",
                                    "21q", "22q"))
            arm.labels = setNames(arm$arm, arm$arm)
            arm.labels[names(small.arms)] = small.arms

            message("Plot Fragmentation profiles (ratio, coverage, combined).")
            ## Plot fragment ratio profile
            f1 =
              ggplot(df.fr3, aes(x=bin, y=ratio_centered, group=sample_id, color="red")) +
              geom_line(size=0.75, alpha=0.75)
            f1 = f1 + geom_line(data=healthy_median, aes(x=bin, y=median_ratio_centered), size=0.75, alpha=0.5, color="black")
            f1 = f1 + labs(x="", y="Fragmentation profile\n", color="")
            f1 = f1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
            f1 = f1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
            f1 = f1 + mytheme
            ggsave(file.path(outdir, paste0(id, "_plot_fragment.pdf")), f1, width=15, height=3, units="in")

            ## Plot short fragment coverage profile
            c1 =
              ggplot(df.fr3, aes(x=bin, y=coverage_centered, group=sample_id, color="red")) +
              geom_line(size=0.75, alpha=0.75)
            c1 = c1 + geom_line(data=healthy_median, aes(x=bin, y=median_coverage_centered), size=0.75, alpha=0.5, color="black")
            c1 = c1 + labs(x="", y="Coverage profile\n", color="")
            c1 = c1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
            c1 = c1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
            c1 = c1 + mytheme
            ggsave(file.path(outdir, paste0(id, "_plot_coverage.pdf")), c1, width=15, height=3, units="in")

            ## Plot combined profile
            b1 =
              ggplot(df.fr3, aes(x=bin, y=combined_centered, group=sample_id, color="red")) +
              geom_line(size=0.75, alpha=0.75)
            b1 = b1 + geom_line(data=healthy_median, aes(x=bin, y=median_combined_centered), size=0.75, alpha=0.5, color="black")
            b1 = b1 + labs(x="", y="Combined profile\n", color="")
            b1 = b1 + facet_grid(~arm, switch="x",space="free_x", scales="free_x", labeller=labeller(arm=arm.labels))
            b1 = b1 + coord_cartesian(xlim = NULL, ylim=c(-.10,.12), expand = TRUE)
            b1 = b1 + mytheme
            ggsave(file.path(outdir, paste0(id, "_plot_combined.pdf")), b1, width=15, height=3, units="in")
          })
