#' @title The fragSet Class
#'
#' @description The fragSet object is a storage container used to store
#'     all intermediate and final results of fragmentomics analysis
#'
#' @slot sampleID Character of sample ID
#' @slot refGenome Character of reference genome used (default is BSgenome.Hsapiens.UCSC.hg38)
#' @slot inputFileType Character of input file type ("bam" or "bedpe")
#' @slot frags List to store non-mitochondrial fragments GRanges object
#' @slot galp_len Numeric of the length of genomic alignment pair vector (i.e. total number of fragments sequenced)
#' @slot mito List to store mitochondrial fragments GRanges object, and mitochondrial fragment report
#' @slot bins List to store 100kb bins and 5Mb bins results
#' @slot results List to store 5Mb distance to healthy median, 5Mb correlation to healthy median, fragmentomics profile plots
#' @slot QC List to store high coverage regions, coverage statistics and GC metrics
#'
#' @exportClass fragSet
#'
setClass(Class = "fragSet",
         slots =  list(sampleID="character",
                       refGenome="character",
                       inputFileType="character",
                       frags="list", ## stores frags_Granges
                       galp_len="numeric",
                       mito="list", ## stores mito_Granges, mitoReport
                       bins="list", ## stores 100kb bins, 5Mb bins
                       results="list", ## stores 5Mb dist, summary, fragment/coverage/combined.pdf
                       QC="list") ## stores high_cov_regions, cov_stats, GC_metrics
)

#' @title Initiate fragSet object
#'
#' @description createFragSet() will read in BAM or BEDPE input file for a sample,
#'     filter reads based on SAM FLAG, fragment width, and mapQ threshold. Then
#'     stores all fragments in a GRanges object with GC content information from
#'     the reference genome specified. This frags GRanges object is stored in the
#'     "frags" slot. Fragments mapping to mitochondrial genome is stored in a
#'     separate GRanges object and stored in the "mito" slot.
#'
#' @param sampleID Character of sample ID
#' @param filePath Character of full path to input file
#' @param fileType Character of either "bam" or "bedpe"
#' @param refGenome Character of reference genome used (default is BSgenome.Hsapiens.UCSC.hg38)
#' @param chr.select Vector of chromosome names (default is paste0("chr",1:22))
#' @param width_min Numeric of minimum fragment length (default is 90)
#' @param width_max Numeric of maximum fragment length (default is 220)
#' @param mapQthreshold Numeric of mapping quality filter greater or equal to value (default is 30)
#' @param keepDuplicates Boolean of whether to keep duplicates reads (default is FALSE - remove duplicate reads with FLAG 0x400 set)
#' @param keepSecondaryAlignments Boolean of whether to keep secondary reads (default is FALSE - remove secondary reads with FLAG 0x100 set)
#'
#' @return Returns a fragSet object
#'
#' @example
#' fset = createFragSet(sampleID = "some_sample",
#'     filePath = "/full/path/to/a_sample.bam",
#'     fileType = "bam",
#'     refGenome="BSgenome.Hsapiens.UCSC.hg38",
#'     chr.select=paste0("chr",1:22),
#'     width_min=90, width_max=220, mapQthreshold=30,
#'     keepDuplicates=FALSE, keepSecondaryAlignments=FALSE)
#'
#' @import Rsamtools
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import BSgenome.Hsapiens.UCSC.hg38
#' @importFrom methods new
#' @importFrom biovizBase GCcontent
#' @importFrom GenomeInfoDb keepSeqlevels
#'
#' @export
#'
#' @include BEDPE_read_and_filter_forFrag.R
#'
createFragSet =
  function(sampleID, filePath, fileType,
           refGenome="BSgenome.Hsapiens.UCSC.hg38",
           chr.select=paste0("chr",1:22),
           width_min=90, width_max=220, mapQthreshold=30,
           keepDuplicates=FALSE, keepSecondaryAlignments=FALSE){
    if (fileType == "bam"){
      indexed.bam = gsub("$", ".bai", filePath)
      if (!file.exists(indexed.bam)) {
        indexBam(filePath)
      }
      param = ScanBamParam(flag = scanBamFlag(isPaired = TRUE,
                                              isProperPair = TRUE,
                                              isDuplicate = keepDuplicates,
                                              isSecondaryAlignment = keepSecondaryAlignments,
                                              isUnmappedQuery = FALSE),
                           mapqFilter = mapQthreshold)
      message("Reading in BAM as genome aligned pairs...")
      galp = readGAlignmentPairs(filePath, param = param)
      rm(param, indexed.bam)

      message("Generating genome fragments GRanges object.")
      frags = granges(GenomeInfoDb::keepSeqlevels(galp, chr.select, pruning.mode="coarse"),
                      on.discordant.seqnames="drop")
      message("Generating mitochondrial GRanges object.")
      mito = granges(keepSeqlevels(galp, paste0("chrM"), pruning.mode="coarse"),
                     on.discordant.seqnames="drop")
      galp_len = length(galp)
      rm(galp)

      w.all = width(frags)
      frags = frags[which(w.all >= width_min & w.all <= width_max)]
      rm(w.all)

      message("Attaching reference GCcontent to fragment object...")
      hsapiens = getBSgenome(refGenome)
      gcs = GCcontent(hsapiens, unstrand(frags))
      frags$gc = gcs
      rm(gcs)

      fs = new('fragSet',
               sampleID=sampleID,
               refGenome=refGenome,
               inputFileType=fileType,
               frags=list("frags_Granges"=frags),
               galp_len=galp_len,
               mito=list("mito_Granges"=mito))

    } else if (fileType == "bedpe"){
      message("Reading BEDPE...")
      bedpe_raw = readBEDPE(filePath)

      message("Filtering BEDPE...")
      bedpe_filtered =
        filtBEDPE(bedpe = bedpe_raw,
                  mapQ = mapQthreshold,
                  isPaired = TRUE,
                  isProperPair = TRUE,
                  isUnmappedQuery = FALSE,
                  isSecondaryAlignment = keepSecondaryAlignments,
                  isDuplicate = keepDuplicates)

      bedpe4fragm = bedpe4frag(bedpe_filtered)
      names(bedpe4fragm) = c("chr", "start", "end", "strand")

      bedpe4fragm.chr1to22 =
        bedpe4fragm[bedpe4fragm$chr %in% chr.select,]

      bedpe4fragm.mito =
        bedpe4fragm[bedpe4fragm$chr == "chrM",]

      message("Generating genome fragments GRanges object.")
      frags = makeGRangesFromDataFrame(df = bedpe4fragm.chr1to22, keep.extra.columns = F)
      message("Generating mitochondrial GRanges object.")
      mito = makeGRangesFromDataFrame(df = bedpe4fragm.mito, keep.extra.columns = F)

      bedpe_len = nrow(bedpe4fragm)
      rm(bedpe4fragm)

      w.all = width(frags)
      frags = frags[which(w.all >= width_min & w.all <= width_max)]
      rm(w.all)

      message("Attaching reference GCcontent to fragment object...")
      hsapiens = getBSgenome(refGenome)
      gcs = GCcontent(hsapiens, unstrand(frags))
      frags$gc = gcs
      rm(gcs)

      fs = new('fragSet',
               sampleID=sampleID,
               refGenome=refGenome,
               inputFileType=fileType,
               frags=list("frags_Granges"=frags),
               galp_len=bedpe_len,
               mito=list("mito_Granges"=mito))
    }
  }
