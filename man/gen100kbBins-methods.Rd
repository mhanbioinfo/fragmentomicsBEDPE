% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frag.gen100kbBins.R
\name{gen100kbBins}
\alias{gen100kbBins}
\alias{gen100kbBins,fragSet-method}
\alias{gc.correct}
\alias{gc.pred}
\title{Generate 100kb bin fragment coverage dataframe}
\usage{
gen100kbBins(
  object,
  filters,
  gaps,
  VNTRs,
  tiles,
  refGenome = "BSgenome.Hsapiens.UCSC.hg38",
  chr.select = paste0("chr", 1:22)
)

\S4method{gen100kbBins}{fragSet}(
  object,
  filters,
  gaps,
  VNTRs,
  tiles,
  refGenome = "BSgenome.Hsapiens.UCSC.hg38",
  chr.select = paste0("chr", 1:22)
)

gc.correct(coverage, bias)

gc.pred(coverage, bias)
}
\arguments{
\item{object}{fragSet object}

\item{filters}{Character of path to filters rda file
(ENCFF356LFX.bed blacklist regions (910 ranges) to filter out)}

\item{gaps}{Character of path to gaps rda file
(centromere, telomere and contig gaps (651 ranges) to filter out)}

\item{VNTRs}{Character of path to VNTRs rda file
(variable number tandem repeats to (382 ranges) to filter out)}

\item{tiles}{Character of path to tiles bed file
(30894 100kb bins template used to put fragments in for coverage calculation)}

\item{refGenome}{Character of reference genome used
(default is "BSgenome.Hsapiens.UCSC.hg38")}

\item{chr.select}{Vector for selecting chromosomes to include in analysis
(default is paste0("chr",1:22))}
}
\value{
Returns a fragSet object with 100kb bin fragment coverage dataframe stored
in "bins" slot, access with get100kbBins() function
}
\description{
The 100kb bin fragment coverage dataframe including the following columns:
"seqnames", "start", "end", "width", "strand", "Seqlength",
"arm", "C.G", "short", "long", "ratio", "nfrags", "coverage",
"short.corrected", "long.corrected", "nfrags.corrected", "ratio.corrected",
"coverage.corrected", "combined", "short.predicted", "long.predicted",
"nfrags.predicted", "ratio.predicted", "coverage.predicted",
"mode", "mean", "median", "quantile.25", "quantile.75", "frag.gc",
"id", and number of fragments with width of 90 up to 220
}
