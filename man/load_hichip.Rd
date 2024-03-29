% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{load_hichip}
\alias{load_hichip}
\title{Load in the hichip matrix as a a sparseMatrix
type specifies what type the file is. Currently only supports the hic-pro
.matrix format. Should be made from allValidPairs.}
\usage{
load_hichip(filename, coords, chrom = NULL, type = "hic-pro", diagonal = 3)
}
\arguments{
\item{filename}{specifies the name of the file to load in}

\item{coords}{specifies the total dimensions of the matrix, created from bins
using load_coords.}

\item{chrom}{specifies which chromosome to subset to based on coords. Default
is NULL so it will return the whole thing.}

\item{type}{specifies the file type}

\item{diagonal}{specifies which entries to 0 out from the diagonal.}
}
\description{
Load in the hichip matrix as a a sparseMatrix
type specifies what type the file is. Currently only supports the hic-pro
.matrix format. Should be made from allValidPairs.
}
