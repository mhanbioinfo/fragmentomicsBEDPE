#' @title Generate 5Mb bin fragment coverage dataframe from 100kb bin dataframe
#'
#' @description
#' The 5Mb bin fragment coverage dataframe including the following columns:
#'     "sample_id", "seqnames", "arm", "start", "end", "gc", "frag_gc",
#'     "short", "long", "nfrags", "ratio", "short_corrected", "long_corrected",
#'     "nfrags_corrected", "ratio_corrected", "ratio_centered", "coverage",
#'     "coverage_corrected", "coverage_centered", "combined", "combined_centered",
#'     "short_var", "long_var", "nfrags_var", "mode_size", "mean_size",
#'     "median_size", "q25_size", "q75_size", "binsize", "bin", "combine"
#'
#' @param object fragSet object
#'
#' @return Returns a fragSet object with 5Mb bin fragment coverage dataframe stored
#'     in "bins" slot, access with get5MbBins() function
#'
#' @example fset = gen5MbBins(fset)
#'
#' @import tidyverse
#' @rdname gen5MbBins-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('gen5MbBins',
           function(object) standardGeneric('gen5MbBins'))
#' @rdname gen5MbBins-methods
setMethod('gen5MbBins', 'fragSet',
          function(object){

            tib.list = dplyr::select(object@bins$kb100bins, -matches("X"))

            df.fr2 = tib.list
            rm(tib.list)
            armlevels = c("1p","1q","2p","2q","3p","3q","4p","4q","5p","5q","6p","6q",
                          "7p","7q","8p","8q", "9p", "9q","10p","10q","11p","11q","12p",
                          "12q","13q","14q","15q","16p","16q","17p","17q","18p","18q",
                          "19p", "19q","20p","20q","21q","22q")
            df.fr2$arm = factor(df.fr2$arm, levels=armlevels)

            message("Combine adjacent 100kb bins to form 5mb bins.")
            df.fr2 =
              df.fr2 %>%
              group_by(arm) %>%
              dplyr::mutate(combine = ifelse(grepl("p", arm),
                                             ceiling((1:length(arm))/50),
                                             ceiling((1:length(arm))/50)))
            df.fr3 =
              df.fr2 %>%
              group_by(id, seqnames, arm, combine) %>%
              dplyr::summarize(short2=sum(short, na.rm=TRUE),
                               long2=sum(long, na.rm=TRUE),
                               short.corrected2=sum(short.corrected, na.rm=TRUE),
                               long.corrected2=sum(long.corrected, na.rm=TRUE),
                               gc=mean(C.G, na.rm=TRUE),
                               frag.gc2=mean(frag.gc, na.rm=TRUE),
                               nfrags2=sum(nfrags, na.rm=TRUE),
                               nfrags.corrected2=sum(nfrags.corrected, na.rm=TRUE),
                               short.var=var(short.corrected, na.rm=TRUE),
                               long.var=var(long.corrected, na.rm=TRUE),
                               nfrags.var=var(nfrags.corrected, na.rm=TRUE),
                               mode_size=unique(mode, na.rm=TRUE),
                               mean_size=unique(mean, na.rm=TRUE),
                               median_size=unique(median, na.rm=TRUE),
                               q25_size=unique(quantile.25, na.rm=TRUE),
                               q75_size=unique(quantile.75, na.rm=TRUE),
                               start=start[1],
                               end=rev(end)[1],
                               binsize=n())

            df.fr3$ratio2 = df.fr3$short2/df.fr3$long2
            df.fr3$ratio.corrected2 = df.fr3$short.corrected2/df.fr3$long.corrected2
            df.fr3$coverage2 = df.fr3$short2/sum(df.fr3$nfrags2)
            df.fr3$coverage.corrected2 = df.fr3$short.corrected2/sum(df.fr3$nfrags.corrected2)
            df.fr3$combined2 = df.fr3$ratio.corrected2*df.fr3$coverage.corrected2

            ## Assign bins
            df.fr3 = df.fr3 %>% dplyr::filter(binsize==50)
            df.fr3 = df.fr3 %>% group_by(id) %>% dplyr::mutate(bin = 1:length(id))

            message("Z-score transform fragment ratio, coverage and combined score.")
            df.fr3$ratio.centered = ((df.fr3$ratio.corrected2 - mean(df.fr3$ratio.corrected2))/sd(df.fr3$ratio.corrected2))*0.01
            df.fr3$coverage.centered = ((df.fr3$coverage.corrected2 - mean(df.fr3$coverage.corrected2))/sd(df.fr3$coverage.corrected2))*0.01
            df.fr3$combined.centered = ((df.fr3$combined2 - mean(df.fr3$combined2))/sd(df.fr3$combined2))*0.01

            ### Rename and reorder dataframe
            df.fr3 =
              df.fr3 %>%
              dplyr::rename(
                sample_id = id,
                short = short2,
                long = long2,
                short_corrected = short.corrected2,
                long_corrected = long.corrected2,
                ratio = ratio2,
                ratio_corrected = ratio.corrected2,
                nfrags = nfrags2,
                nfrags_corrected = nfrags.corrected2,
                coverage_corrected = coverage.corrected2,
                short_var = short.var,
                long_var = long.var,
                nfrags_var = nfrags.var,
                ratio_centered = ratio.centered,
                coverage_centered = coverage.centered,
                combined = combined2,
                combined_centered = combined.centered,
                frag_gc = frag.gc2,
                coverage = coverage2
              )

            df.fr3 =
              df.fr3 %>%
              dplyr::relocate(sample_id, seqnames, arm, start, end, gc, frag_gc, short, long, nfrags, ratio,
                              short_corrected, long_corrected, nfrags_corrected, ratio_corrected, ratio_centered,
                              coverage, coverage_corrected, coverage_centered, combined, combined_centered,
                              short_var, long_var, nfrags_var, mode_size,mean_size, median_size, q25_size, q75_size,
                              binsize, bin)

            object@bins$mb5bins = df.fr3
            rm(df.fr3)

            return(object)
          })
