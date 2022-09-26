#' @title Compare sample fragmentation profile to healthy median
#'
#' @description
#' Calculate correlation and distance of sample fragmentation profile to
#'     healthy median (provided with package or generated externally).
#'
#' @param object fragSet object
#' @param healthy Character of path to healthy median rda file
#'
#' @return Returns a fragSet object with mb5dist and md5summary dataframes
#'     stored in results slot, access with get5MbDist() and get5MbSummary()
#'     respectively
#'
#' @example fset = compareToHealthy(fset, healthy=healthy_path)
#'
#' @import tidyverse
#' @rdname compareToHealthy-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('compareToHealthy',
           function(object, healthy) standardGeneric('compareToHealthy'))
#' @rdname compareToHealthy-methods
setMethod('compareToHealthy', 'fragSet',
          function(object, healthy){

            load(healthy)

            sd_dist = function(sample, mean, sd) {
              sds = (sample - mean)/sd
            }

            df.fr3 = object@bins$mb5bins

            message("Generate correlations between sample and healthy control.")
            correlations = df.fr3 %>% ungroup() %>% group_by(sample_id) %>%
              dplyr::summarize(ratio_corr=cor(ratio, healthy_median$median_ratio, method="pearson", use="complete.obs"),
                               ratio_corrected_corr=cor(ratio_corrected, healthy_median$median_ratio_corrected, method="pearson", use="complete.obs"),
                               ratio_centered_corr=cor(ratio_centered, healthy_median$median_ratio_centered, method="pearson", use="complete.obs"),
                               coverage_corr=cor(coverage, healthy_median$median_coverage, method="pearson", use="complete.obs"),
                               coverage_corrected_corr=cor(coverage_corrected, healthy_median$median_coverage_corrected, method="pearson", use="complete.obs"),
                               coverage_centered_corr=cor(coverage_centered, healthy_median$median_coverage_centered, method="pearson", use="complete.obs"),
                               combined_corr=cor(combined, healthy_median$median_combined, method="pearson", use="complete.obs"),
                               combined_centered_corr=cor(combined_centered, healthy_median$median_combined_centered, method="pearson", use="complete.obs"),
                               nfrags = sum(nfrags),
                               mode_size=unique(mode_size),
                               mean_size=unique(mean_size),
                               median_size=unique(median_size),
                               q25_size=unique(q25_size),
                               q75_size=unique(q75_size),
                               hqbases_analyzed=100*sum(nfrags)*2,
                               depth=hqbases_analyzed/(504*5e6)
              )

            message("Generate distance between sample and healthy control.")
            distance = df.fr3 %>%
              dplyr::summarize(seqnames=seqnames,
                               arm=arm,
                               start=start,
                               end=end,
                               gc=gc,
                               ratio_dist=sd_dist(ratio, healthy_median$median_ratio, healthy_median$sd_ratio),
                               ratio_corrected_dist=sd_dist(ratio_corrected, healthy_median$median_ratio_corrected, healthy_median$sd_ratio_corrected),
                               ratio_centered_dist=sd_dist(ratio_centered, healthy_median$median_ratio_centered, healthy_median$sd_ratio_centered),
                               coverage_dist=sd_dist(coverage, healthy_median$median_coverage, healthy_median$sd_coverage),
                               coverage_corrected_dist=sd_dist(coverage_corrected, healthy_median$median_coverage_corrected, healthy_median$sd_coverage_corrected),
                               coverage_centered_dist=sd_dist(coverage_centered, healthy_median$median_coverage_centered, healthy_median$sd_coverage_centered),
                               combined_dist=sd_dist(combined, healthy_median$median_combined, healthy_median$sd_combined),
                               combined_centered_dist=sd_dist(combined_centered, healthy_median$median_combined_centered, healthy_median$sd_combined_centered)
              )

            message("Generate summary of correlation and distance between sample and healthy control.")
            summary_sd = distance %>% ungroup() %>% group_by(sample_id) %>%
              dplyr::summarize(ratio_dist=mean(ratio_dist),
                               ratio_corrected_dist=mean(ratio_corrected_dist),
                               ratio_centered_dist=mean(ratio_centered_dist),
                               coverage_dist=mean(coverage_dist),
                               coverage_corrected_dist=mean(coverage_corrected_dist),
                               coverage_centered_dist=mean(coverage_centered_dist),
                               combined_dist=mean(combined_dist),
                               combined_centered_dist=mean(combined_centered_dist)
              )

            summary_df =
              data.frame(metrics = c(names(unlist(correlations)), names(unlist(summary_sd))),
                         value = c(unlist(correlations), unlist(summary_sd))) %>%
              dplyr::distinct(metrics, .keep_all=T)

            object@results$mb5dist = distance
            object@results$md5summary = summary_df

            rm(correlations, distance, summary_df, sd_dist)

            return(object)
          })
