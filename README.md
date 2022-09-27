# fragmentomicsBEDPE

## Description

The fragmentomicsBEDPE R package peforms fragmentomics analysis on whole genome sequencing data. Fragmentomics analysis assess 2 main metrics: the short to long fragment ratio (hereafter as “ratio”), short fragment coverage (hereafter as “coverage”), as well as combination of ratio*coverage (hereafter as “combined”) in 5Mb bins throughout the entire genome.

The package then compares a sample’s ratio, coverage and combined metrics against healthy median, and generate 2 additional metrics: correlation to healthy median, and distance to healthy median.

This package packages scripts from Derek Wong’s fragmentomics pipeline (https://github.com/derekwong90/fragmentomics/), which was based on DELPHI (https://github.com/cancer-genomics/delfi_scripts).

## Installation

```
require(devtools)
require(tidyverse)
devtools::install_github("mhanbioinfo/fragmentomicsBEDPE")

library(fragmentomicsBEDPE)
```
Please follow vignette for example run (https://github.com/mhanbioinfo/fragmentomicsBEDPE/blob/main/vignettes/fragmentomicsBEDPE_vignette.Rmd)
