\name{pcaReduce-package}
\alias{pcaReduce-package}
\alias{pcaReduce}
\docType{package}
\title{pcaReduce}
\description{Hierarchical representation  of Single Cell RNA-seq data via dimensionality reduction and clustering}
\details{
\tabular{ll}{
Package: \tab pcaReduce\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-04-28\cr
License: \tab GPL (>= 2)\cr
}
Developed to be used on Single Cell RNA-seq data (but not limited to SC data). The data should be presented in a matrix format
}
\author{Justina Zurauskiene and Christopher Yau

Maintainer: Justina Zurauskiene <justina_zurauskiene@gmail.com>
}
\references{
Justina Zurauskiene, Christopher Yau (2015). ``pcaReduce: Hierarchical Clustering of Single Cell Transcriptional Profiles". Journal XXXX, xx(xx), xx-xx. URL: http://xxxxx.

}
~~ Optionally other standard keywords, one per line, from file ~~
~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
# This is example how to use it:
# d_test <- read.delim("yourDataset.csv", sep=",", header=T)
# Names <- read.table("Annotations_artificial_data.csv", sep=",", header = TRUE)

# rownames(d_test) <- Names[,1]

# repear 5 times, project initialy to 10 dimensions, "M" merge based on max:
# L <- PCAreduce(D_t, 5, 10, "M")

# repear 5 times, project initialy to 10 dimensions, "S" merge based on sampling:
# L <- PCAreduce(D_t, 5, 10, "S")
}
