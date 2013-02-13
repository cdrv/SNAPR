\name{LDSearch}
\alias{LDSearch}
\title{Search for SNPs in LD with a set of SNPs}
\usage{
  LDSearch(SNPs, dataset = "onekgpilot", panel = "CEU",
    RSquaredLimit = 0.8, distanceLimit = 500)
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
  \item{\code{YRI}} \item{\code{JPT+CHB}} } If you are
  working with \code{hapmap3r2}, you can choose the
  additional panels: \itemize{ \item{\code{ASW}}
  \item{\code{CHD}} \item{\code{GIH}} \item{\code{LWK}}
  \item{\code{MEK}} \item{\code{MKK}} \item{\code{TSI}}
  \item{\code{CEU+TSI}} \item{\code{JPT+CHB+CHD}} }}

  \item{RSquaredLimit}{The R Squared limit to specify as a
  filter for returned SNPs; that is, only SNP pairs with
  R-squared greater than \code{RSquaredLimit} will be
  returned.}

  \item{distanceLimit}{The distance (in kilobases) upstream
  and downstream to search for SNPs in LD with each set of
  SNPs.}
}
\value{
  A list of data frames, one for each SNP queried,
  containing information about the SNPs found to be in LD
  with that SNP.
}
\description{
  This function queries SNAP for SNPs in high linkage
  disequilibrium with a set of SNPs.
}
\details{
  For more details, please see
  \url{http://www.broadinstitute.org/mpg/snap/ldsearch.php}.

  Information on the HapMap populations:
  \url{http://ccr.coriell.org/Sections/Collections/NHGRI/hapmap.aspx?PgId=266&coll=HG}

  Information on the 1000 Genomes populations:
  \url{http://www.1000genomes.org/category/frequently-asked-questions/population}
}
\examples{
\dontrun{LDSearch("rs429358")}
}

