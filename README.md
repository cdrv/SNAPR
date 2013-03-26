SNAPR
=====
  
This package provides an R interface to SNP Annotation and Proxy Search, as
hosted by the Broad Institute [here](http://www.broadinstitute.org/mpg/snap/ldsearch.php).
In particular, it exposes an R interface to the `ldsearch.php` application.

Install me with `devtools`:

    install_github("SNAPR", "cdrv")
    
After installation, you can try a simple query with

    library(SNAPR)
    LDSearch("rs420358")
    
