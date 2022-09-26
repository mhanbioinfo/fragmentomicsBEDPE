#' @title Getter functions related to fset object
#'
#' @description
#' Getter functions to access slots of fset object
#'
#' @param object fragSet object
#'
#' @return various
#'
#' @rdname fsetGetter-methods
#'
#' @include frag.fragSetClass_createFragSet.R
NULL

#' @rdname fsetGetter-methods
#' @export
setMethod('show', signature='fragSet', definition=function(object) {
  message("fragSet")
  message("=======================================")
  message(paste0("Sample ID: ", getSampleID(object), " ."))
})

#' @rdname fsetGetter-methods
#' @export
setGeneric('getSampleID',
           function(object) standardGeneric('getSampleID'))
#' @rdname fsetGetter-methods
setMethod('getSampleID', 'fragSet',
          function(object){
            sampleID=object@sampleID
            return(sampleID)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getGalpLen',
           function(object) standardGeneric('getGalpLen'))
#' @rdname fsetGetter-methods
setMethod('getGalpLen', 'fragSet',
          function(object){
            galp_len=object@galp_len
            return(galp_len)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getMitoGRanges',
           function(object) standardGeneric('getMitoGRanges'))
#' @rdname fsetGetter-methods
setMethod('getMitoGRanges', 'fragSet',
          function(object){
            mito_Granges=object@mito$mito_Granges
            return(mito_Granges)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getFragsGRanges',
           function(object) standardGeneric('getFragsGRanges'))
#' @rdname fsetGetter-methods
setMethod('getFragsGRanges', 'fragSet',
          function(object){
            frags_Granges=object@frags$frags_Granges
            return(frags_Granges)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getMitoReport',
           function(object) standardGeneric('getMitoReport'))
#' @rdname fsetGetter-methods
setMethod('getMitoReport', 'fragSet',
          function(object){
            mt_matrix = object@mito$mitoReport
            return(mt_matrix)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getHighCovRegions',
           function(object) standardGeneric('getHighCovRegions'))
#' @rdname fsetGetter-methods
setMethod('getHighCovRegions', 'fragSet',
          function(object){
            high_cov_regions = object@QC$high_cov_regions
            return(high_cov_regions)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('getCovStats',
           function(object) standardGeneric('getCovStats'))
#' @rdname fsetGetter-methods
setMethod('getCovStats', 'fragSet',
          function(object){
            cov_stats = object@QC$cov_stats
            return(cov_stats)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('get100kbBins',
           function(object) standardGeneric('get100kbBins'))
#' @rdname fsetGetter-methods
setMethod('get100kbBins', 'fragSet',
          function(object){
            tib.list = object@bins$kb100bins
            return(tib.list)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('get5MbBins',
           function(object) standardGeneric('get5MbBins'))
#' @rdname fsetGetter-methods
setMethod('get5MbBins', 'fragSet',
          function(object){
            df.fr3 = object@bins$mb5bins
            return(df.fr3)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('get5MbDist',
           function(object) standardGeneric('get5MbDist'))
#' @rdname fsetGetter-methods
setMethod('get5MbDist', 'fragSet',
          function(object){
            distance = object@results$mb5dist
            return(distance)
          })

#' @rdname fsetGetter-methods
#' @export
setGeneric('get5MbSummary',
           function(object) standardGeneric('get5MbSummary'))
#' @rdname fsetGetter-methods
setMethod('get5MbSummary', 'fragSet',
          function(object){
            summary_df = object@results$md5summary
            return(summary_df)
          })
