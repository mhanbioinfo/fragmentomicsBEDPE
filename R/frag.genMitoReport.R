#' @title Generate mitochondrial fragment report
#'
#' @description
#' Mitochondrial fragment report includes: mt_short (number of mito short reads),
#'     mt_long (number of mito long reads), mt_nFrag (total number of fragments mapped to mito),
#'     mt_ratio (ratio of number of short/long mito reads), mt_coverage (sum(mt_width)/16569),
#'     mt_fraction (number of mitochondrial fragments / total number of fragments),
#'     mt_median/mean/mode/lower_quartile/upper_quartile
#'     (median/mean/mode/lower and upper quartile of mitochondrial fragment distribution),
#'     and sampleID
#'
#' @param object fragSet object
#'
#' @return Returns a fragSet object with mitochondrial fragment report stored in "mito" slot, access with getMitoReport() function
#'
#' @example fset = genMitoReport(fset)
#'
#' @importFrom stats density
#' @importFrom matrixStats count
#' @rdname genMitoReport-methods
#' @export
#'
#' @include frag.showAndGetterFunctions.R
#'
setGeneric('genMitoReport',
           function(object) standardGeneric('genMitoReport'))
#' @rdname genMitoReport-methods
setMethod('genMitoReport', 'fragSet',
          function(object){
            sampleID=getSampleID(object)
            mito_Granges=getMitoGRanges(object)
            galp_len=getGalpLen(object)

            ## Extract Mitochondrial reads
            mt_nFrag = length(mito_Granges)
            mt_width = width(mito_Granges)
            rm(mito_Granges)

            ## Calculate mitochondrial stats
            mt_median = median(mt_width)
            mt_mean = mean(mt_width)
            dens = stats::density(mt_width)
            mt_mode = dens$x[which.max(dens$y)]
            mt_lower_quartile = as.integer(quantile(mt_width, 0.25))
            mt_upper_quartile = as.integer(quantile(mt_width, 0.75))
            mt_short = matrixStats::count(mt_width <= 150)
            mt_long = matrixStats::count(mt_width > 150)
            mt_ratio = mt_short/mt_long
            mt_coverage = sum(mt_width)/16569
            mt_fraction = mt_nFrag/galp_len

            ## Generate mitochondrial report
            mt_matrix = data.frame(mt_short, mt_long, mt_nFrag, mt_ratio, mt_coverage, mt_fraction,
                                   mt_median, mt_mean, mt_mode, mt_lower_quartile, mt_upper_quartile, sampleID)
            object@mito$mitoReport = mt_matrix
            rm(mt_median, mt_mean, dens, mt_mode, mt_lower_quartile, mt_upper_quartile, mt_short, mt_long,
               mt_ratio, mt_coverage, mt_fraction, mt_nFrag, mt_width, mt_matrix)

            return(object)
          })
