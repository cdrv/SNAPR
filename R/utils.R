#' Swap Elements in a Vector
#' 
#' @param vec A character vector, or vector coercable to character.
#' @param from A vector of elements to map from.
#' @param to A vector of elements to map to.
#' @param ... Optional arguments passed to \code{\link{match}}.
swap <- function( vec, from, to=names(from), ... ) {
  tmp <- to[ match( vec, from, ... ) ]
  tmp[ is.na(tmp) ] <- vec[ is.na(tmp) ]
  return( tmp )
}

#' Flip Genotypes
#' 
#' Given a set of genotypes, flip them.
#' @param SNP A vector of genotypes for a particular locus.
#' @param sep The separator between each allele.
#' @param outSep The output separator to use.
flip <- function( SNP, sep="", outSep=sep ) {
  
  base1 <- c("A","C","G","T")
  base2 <- c("T","G","C","A")
  
  from <- do.call( 
    function(...) { paste(..., sep=sep) }, 
    expand.grid( 
      base1, base1
    ) )
  
  to <- do.call( function(...) { paste(..., sep=outSep) }, expand.grid( 
    base2, base2
  ) )
  
  return( swap( SNP, from, to ) )
  
}

#' Query NCBI for info on a set of SNPs
#' 
#' @param SNPs A vector of SNPs (rs numbers).
get_snp_info <- function(SNPs) {
  
  ## Query the unmatched SNPs
  query <- "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&mode=xml&id="
  query <- paste0( query, paste( SNPs, collapse="," ) )
  
  xml <- getURL( query )
  xml_parsed <- xmlInternalTreeParse( xml )
  xml_list <- xmlToList( xml_parsed )
  
  ## The length of the xml list is equal to (num SNPs) + 1, as the last entry
  ## is the schema location. So we iterate through items 1 to (n-1).'
  
  out <- NULL
  
  for( i in 1:(length(xml_list)-1) ) {
    
    myList <- xml_list[[i]]
    my_chr <- myList$Assembly$Component$.attrs["chromosome"]
    my_snp <- paste0("rs", myList$.attrs["rsId"])
    
    my_gene <- tryCatch( 
      myList$Assembly$Component$MapLoc$FxnSet['symbol'],
      error = function(e) { return(NA) }
    )
    
    ## by default, we will take data from the forward strand
    ## if we have grabbed something from the reverse, we flip it
    
    tmp <- c( myList$Ss$Sequence$Observed, myList$Ss$.attrs["orient"] )
    if( tmp[2] != "forward" ) {
      tmp[1] <- flip( tmp[1], sep="/" )
    }
    
    alleles <- unlist( strsplit( tmp[1], "/" ) )
    
    ## check which of the two alleles grabbed is actually the minor allele
    ## we might have to 'flip' if there is no match
    maf_allele <- myList$Frequency['allele']
    if( all(maf_allele != alleles) ) {
      maf_allele <- SNAPR:::swap( maf_allele, c("A", "C", "G", "T"), c("T", "G", "C", "A") )
    }
    
    my_minor <- alleles[ maf_allele == alleles ]
    my_major <- alleles[ maf_allele != alleles ]
    my_freq <- myList$Frequency["freq"]
    my_pos <- tryCatch(
      myList$Assembly$Component$MapLoc$.attrs["physMapInt"],
      error=function(e) { 
        myList$Assembly$Component$MapLoc["physMapInt"]
      }
    )
    
    dat <- data.frame( 
      chr = as.integer(my_chr),
      marker = my_snp,
      gene = my_gene,
      major = my_major,
      minor = my_minor,
      maf = as.numeric(my_freq),
      pos = as.integer(my_pos),
      stringsAsFactors=FALSE
    )
    
    out <- rbind( out, dat )      
    
  }
  
  names(out) <- c("chr", "marker", "gene_ncbi", "major_fwd", "minor_fwd", "maf_ncbi", "pos")
  return(out)
  
}

#' Split a Vector of Strings Following a Regular Structure
#' 
#' This function takes a vector of strings following a regular 
#' structure, and converts that vector into a \code{data.frame}, split
#' on that delimiter. A nice wrapper to \code{\link{strsplit}}, essentially 
#' - the primary bonus is the automatic coersion to a \code{data.frame}.
#' @param x a vector of strings.
#' @param sep the delimiter / \code{\link{regex}} you wish to split your strings on.
#' @param fixed logical. If \code{TRUE}, we match \code{sep} exactly; 
#' otherwise, we use regular expressions. Has priority over \code{perl}.
#' @param perl logical. Should perl-compatible regexps be used?
#' @param useBytes logical. If \code{TRUE}, matching is done byte-by-byte rather than
#' character-by-character.
#' @param names optional: a vector of names to pass to the returned \code{data.frame}.
#' @seealso \code{\link{strsplit}}
str_split <- function(x, sep, fixed=FALSE, perl=TRUE, useBytes=FALSE, names=NULL) {
  
  x <- as.character(x)
  
  if( fixed ) {
    perl <- FALSE
  }
  
  tmp <- strsplit( x, sep, fixed=fixed, perl=perl, useBytes=useBytes )
  if( length( unique( unlist( lapply( tmp, length ) ) ) ) > 1 ) {
    stop("non-equal lengths for each entry of x post-splitting")
  }
  tmp <- unlist( tmp )
  tmp <- as.data.frame( 
    matrix( tmp, ncol = (length(tmp) / length(x)), byrow=T ),
    stringsAsFactors=FALSE, optional=TRUE 
  )
  
  if( !is.null(names) ) {
    names(tmp) <- names
  } else {
    names(tmp) <- paste( "V", 1:ncol(tmp), sep="" )
  }
  
  return(tmp)
}
