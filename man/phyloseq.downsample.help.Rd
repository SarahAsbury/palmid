% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phyloseq.downsample.help.R
\name{phyloseq.downsample.help}
\alias{phyloseq.downsample.help}
\title{phyloseq.downsample.help}
\usage{
phyloseq.downsample.help(dat, disaggregated_otu, downsample)
}
\arguments{
\item{dat}{phyloseq object derived from phyloseq.palmid}

\item{disaggregated_otu}{OTU table converted as list with taxa repeated n times. Produced by phyloseq.downsample.}

\item{downsample}{Numeric. Downsample to library depth x.}
}
\value{
downsampled otu table (matrix)
}
\description{
Helper function for phyloseq.downsample. Performs the downsample for a single dataset
}
\seealso{
[phyloseq.palmid, phyloseq.downsample]
}
