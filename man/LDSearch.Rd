\name{LDSearch}
\alias{LDSearch}
\title{Search for SNPs in LD with a set of SNPs}
\usage{
  LDSearch(SNPs, dataset = "onekgpilot", panel = "CEU",
    RSquaredLimit = 0.6, distanceLimit = 10)
}
\arguments{
  \item{SNPs}{A vector of SNPs (rs numbers).}

  \item{dataset}{The dataset to query. Must be one of:
  \itemize{ \item{\code{rel21: }}{HapMap Release 21}
  \item{\code{rel22: }}{HapMap Release 22}
  \item{\code{hapmap3r2: }}{HapMap 3 (release 2)}
  \item{\code{onekgpilot: }}{1000 Genomes Pilot 1} }}

  \item{panel}{The panel to use from the queried data set.
  Must be one of: \itemize{ \item{\code{CEU}}
  \item{\code{YRI}} \item{\code{JPT+CHB}} }}

  \item{RSquaredLimit}{The R Squared limit to specify as a
  filter for returned SNPs; that is, only SNP pairs with
  R-squared greater than \code{RSquaredLimit} will be
  returned.}

  \item{distanceLimit}{The distance (in kilobases) upstream
  and downstream to search for SNPs in LD with each set of
  SNPs.}
}
\value{
  A list of data.frames, one for each SNP queried,
  containing information about the SNPs found to be in LD
  with that SNP.
}
\description{
  This function queries SNAP for SNPs in high linkage
  disequilibrium with a set of SNPs.
}

