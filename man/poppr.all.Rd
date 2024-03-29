% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Index_calculations.r
\name{poppr.all}
\alias{poppr.all}
\title{Process a list of files with poppr}
\usage{
poppr.all(filelist, ...)
}
\arguments{
\item{filelist}{a list of files in the current working directory}

\item{...}{arguments passed on to poppr}
}
\value{
see [poppr()]
}
\description{
poppr.all is a wrapper function that will loop through a list of files from
the working directory, execute [poppr()], and concatenate the
output into one data frame.
}
\examples{
\dontrun{
# Obtain a list of fstat files from a directory.
x <- getfile(multi=TRUE, pattern="^.+?dat$")

# run the analysis on each file.
poppr.all(file.path(x$path, x$files))
}
}
\seealso{
[poppr()], [getfile()]
}
\author{
Zhian N. Kamvar
}
