#' @title Reading and filtering BEDPE file for fragmentomics analysis
#'
#' @description Where readGAlignmentPairs() function from the GenomicAlignments
#'     package reads in BAM file, the 3 functions: readBEDPE(), filtBEDPE() and
#'     bedpe4frag() mimics readGAlignmentPairs() to read in BEDPE file with
#'     identical results.
#'
#' @rdname BEDPE_read_filt
#' @importFrom data.table fread
#' @export
#'
NULL

#' @rdname BEDPE_read_filt
#' @export
readBEDPE = function(bedpe_path){
  BEDPE =
    data.table::fread(file = bedpe_path, sep = '\t', header = F,
                      colClasses = c("character","integer","integer",
                                     "character","integer","integer",
                                     "character","integer","integer",
                                     "character","character","character","character",
                                     "integer","integer","integer","integer","integer","integer")) %>%
    filter(!((bitwAnd(0x100,V14)==0x100 & bitwAnd(0x100,V15)!=0x100) |
               (bitwAnd(0x100,V14)!=0x100 & bitwAnd(0x100,V15)==0x100) |
               (bitwAnd(0x800,V14)==0x800 & bitwAnd(0x800,V15)!=0x800) |
               (bitwAnd(0x800,V14)!=0x800 & bitwAnd(0x800,V15)==0x800))) %>%
    filter(!(bitwAnd(0x4,V14)==0x4 |
               bitwAnd(0x8,V14)==0x8 |
               bitwAnd(0x4,V15)==0x4 |
               bitwAnd(0x8,V15)==0x8))
  ## mandatory removal of frags where one read is secondar or supplementary while the other is not
  ## which readGAlignmentPairs() and Rsamtools in general disregards
  return(BEDPE)
}

#' @rdname BEDPE_read_filt
#' @export
## filtBEDPE functions mimics Rsamtools,
## like Rsamtools, it does not have the full filtering capability of samtools
## i.e. combining flags in AND and OR fashion does not work like samtools
filtBEDPE = function(bedpe,
                     mapQ = NULL,
                     isPaired = NULL,
                     isProperPair = NULL,
                     isUnmappedQuery = FALSE,
                     isSecondaryAlignment = NULL,
                     isDuplicate = NULL,
                     isSupplementaryAlignment = NULL,
                     isTandems = NULL,
                     isR1R2orR2R1 = NULL,
                     isOneSecondaryOtherNot = NULL,
                     hasUnmappedMate = NULL,
                     isMinusStrand = NULL,
                     isMateMinusStrand = NULL,
                     isFirstMateRead = NULL,
                     isSecondMateRead = NULL,
                     isNotPrimaryRead = NULL,
                     isNotPassingQualityControls = NULL) {
  BEDPEfilt = bedpe

  if (!is.null(mapQ)){
    BEDPEfilt =
      BEDPEfilt %>%
      filter(V8 >= mapQ & V9 >= mapQ)
  }
  if (!is.null(isPaired) && isPaired == TRUE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x1,V14)==0x1 & bitwAnd(0x1,V15)==0x1)
  }
  if (!is.null(isProperPair) && isProperPair == TRUE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x2,V14)==0x2)
  } else if (!is.null(isProperPair) && isProperPair == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x2,V14)!=0x2)
  }
  if (!is.null(isUnmappedQuery) && isUnmappedQuery == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter(!(bitwAnd(0x4,V14)==0x4 |
                 bitwAnd(0x8,V14)==0x8 |
                 bitwAnd(0x4,V15)==0x4 |
                 bitwAnd(0x8,V15)==0x8))
  }
  if (!is.null(isSecondaryAlignment) && isSecondaryAlignment == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x100,V14)!=0x100 & bitwAnd(0x100,V15)!=0x100)
  } else if (!is.null(isSecondaryAlignment) && isSecondaryAlignment == TRUE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x100,V14)==0x100 | bitwAnd(0x100,V15)==0x100)
  }
  if (!is.null(isDuplicate) && isDuplicate == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter(bitwAnd(0x400,V14)!=0x400 &
               bitwAnd(0x400,V15)!=0x400)
  }
  if (!is.null(isSupplementaryAlignment) && isSupplementaryAlignment == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter(V1 == V4) %>%
      filter(!((bitwAnd(0x800,V14)==0x800 |
                  bitwAnd(0x800,V15)==0x800) &
                 V1 == V4))
  } else if(!is.null(isSupplementaryAlignment) && isSupplementaryAlignment == TRUE){
    BEDPEfilt = BEDPEfilt %>%
      filter(V1 != V4)
  }

  ## not sure if these are needed
  if (!is.null(isTandems) && isTandems == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter((bitwAnd(0x10,V14)==0x10 &
                bitwAnd(0x20,V15)==0x20) |
               (bitwAnd(0x20,V14)==0x20 &
                  bitwAnd(0x10,V15)==0x10)) %>%
      filter(!((V14==113 & V15==177) |
                 (V14==177 & V15==113)))
  }
  if (!is.null(isR1R2orR2R1) && isR1R2orR2R1 == TRUE){
    BEDPEfilt = BEDPEfilt %>%
      filter((bitwAnd(0x40,V14)==0x40 &
                bitwAnd(0x80,V15)==0x80) |
               (bitwAnd(0x80,V14)==0x80 &
                  bitwAnd(0x40,V15)==0x40))
  }
  if (!is.null(isOneSecondaryOtherNot) && isOneSecondaryOtherNot == FALSE){
    BEDPEfilt = BEDPEfilt %>%
      filter((bitwAnd(0x100,V14)!=0x100 &
                bitwAnd(0x100,V15)!=0x100) |
               (bitwAnd(0x100,V14)==0x100 &
                  bitwAnd(0x100,V15)==0x100))
  }
  return(BEDPEfilt)
}

#' @rdname BEDPE_read_filt
#' @export
bedpe4frag = function(BEDPEfiltd){
  df =
    BEDPEfiltd %>%
    mutate(V2 = V2+1,                       ## convert to 1-base
           V5 = V5+1) %>%
    mutate(V2 = as.integer(V2)) %>%
    mutate(V5 = as.integer(V5)) %>%
    select(V1, V10, V2, V3, V5, V6) %>%
    mutate(start = if_else(V10 == "+", V2, if_else(V10 == "-", V5, as.integer(0))),
           end = if_else(V10 == "+", V6, if_else(V10 == "-", V3, as.integer(0)))) %>%
    select(V1, start, end, V10)
  return(df)
}
