% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2014 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
%               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
%               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
%
\name{GenometriCorrelation}
\alias{GenometriCorrelation}
\alias{GenometricCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
GenometriCorrelation function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
The GenometriCorrelation function compares two interval annotations on a chromosome, set of chromosomes or on an entire genome, and performs various statistical tests to determine whether the two interval sets are independent or are positioned nonrandomly with respect to each other. For a complete description of the tests, refer to the \href{../doc/GenometriCorr.pdf}{the package vignette}.
}
\usage{
	GenometriCorrelation (
		query,reference,
		chromosomes.to.proceed=c(),
		chromosomes.to.include.in.awhole=c(),
		chromosomes.to.exclude.from.awhole=c(),
		add.chr.as.prefix=FALSE,
		awhole.only=FALSE,
		space=NA,
		map.to.half=TRUE,
		showProgressBar=TRUE,
		showTkProgressBar=FALSE,
		chromosomes.length=c(),
		suppress.evaluated.length.warning=FALSE,
		cut.all.over.length=FALSE,
		ecdf.area.permut.number=100,
		mean.distance.permut.number=100,
		jaccard.measure.permut.number=100,
		jaccard.permut.is.rearrangement=FALSE,
		awhole.space.name="awhole",
		keep.distributions=FALSE,
		representing.point.function=mitl,
		query.representing.point.function=representing.point.function,
		reference.representing.point.function=representing.point.function
		)
}

\arguments{
  \item{query}{
The \code{IRanges} or \code{RangedData} or \code{GRanges} object that contains the query interval set coordinates and spaces (generally, chromosomes). The \code{GenometriCorrelation} function tests whether it is positioned independently relative to the \code{reference} interval set.
}
  \item{reference}{
The \code{IRanges} or \code{RangedData} or \code{GRanges} object that contains the reference interval set coordinates and spaces. The \code{GenometriCorrelation} function tests whether \code{query} is positioned independently relative to \code{reference}.
}
  \item{chromosomes.to.proceed}{
This vector of strings contains the names of spaces (chromosomes) to analyze. If both \code{query} and \code{reference} are \code{RangedData} or \code{GRanges} if the parameter is not given (its default is c()), the initial list of spaces to proceed is the intersection of the space lists of query and reference names, if \code{chromosomes.to.proceed} is given, it restricts the list so that only those names from the intersection that are in \code{chromosomes.to.proceed} are analyzed. If one or both of the  \code{query} and \code{reference} is an IRange, only the first element of \code{chromosomes.to.proceed} is used. It actually overrides the \code{space} argument (see below for explanation). If both \code{chromosomes.to.proceed} and \code{space} are given and \code{space != chromosomes.to.proceed[1]}, an error is generated. If after parsing all arguments, less than \code{2} chromosomes are included in the analysis, the overall (awhole) statistics will not be produced.
}
  \item{chromosomes.to.include.in.awhole}{
This vector of strings contains the names of spaces to be included in the overall (awhole) statistics. Its default is c(), meaning that all the analysed genes are included.}
  \item{chromosomes.to.exclude.from.awhole}{
List of chromosomes (spaces) to be excluded from the overall statistics.
}
	\item{add.chr.as.prefix}{The correlation is only performed on chromosomes that have exactly the same name, so by default, a chromosome named chr1 will not be considered the same chromosome as one simply labeled 1. This argument is provided so that if the chromosome names in the \code{query}, \code{reference} and the \code{chromosomes*} parameters, or the names of the chromosomes in \code{chromosome.length}) have no prefix \code{'chr'} (lower-, upper-, or mixed case), this prefix will be added, and all strings matching \code{'chr'} in uppercase or partly uppercase are changed to lowercase. This way, comparisons can be made even if the chromosome names differ by any variation of the word \code{'chr'}. The default is FALSE.
}
  \item{awhole.only}{
If FALSE, all the considered chromosomes statistics and the summary (awhole) information are returned (default). If TRUE, only the summary (awhole) information is returned. Useful for pipelines.
}
  \item{space}{
If both \code{query} and \code{reference} are \code{RangedData}, the parameter is not too useful. It works as \code{chromosomes.to.proceed} of length 1. If one of the two first paramtetesrs is \code{IRanges}, it shows what space from the other one is to be used in the procedure (only one can be used). If both are, it gives the name of the space to be shown in the result.
}
  \item{map.to.half}{
Some of the tests we use are besed on distances between a query point and the closest reference point. If map.to.half is TRUE (default) we look for the closest reference point upstream or downstream, if it is \code{FALSE}, we look only downstream (left). Useful if you are mapping to intervals and want to preserve directionality.
}
  \item{showProgressBar}{
Toggle the text progress bar. Default is \code{TRUE}.
}
  \item{showTkProgressBar}{
Toggle the Tk progress bar. If it is \code{TRUE} but the Tcl/Tk is not loadable (e.g. \code{require('tcltk')} returns \code{FALSE}), the parameter is turned to \code{FALSE}. Default is \code{FALSE}.
}
  \item{chromosomes.length}{
A vector of lengths of the chromosomes to be tested; each chromosome is given a name and a numerical length. The order of the chromosomes does not matter.
}
	\item{suppress.evaluated.length.warning}{
If there is a chromosome that is included in the evaluation but its length is not given in \code{chromosome.length}, the length of the chromosome is calculated as the maximal position mentioned in the data. If this parameter is \code{FALSE}, this will generate a warning, which is suppressed when the parameter is \code{TRUE}. The default is \code{FALSE}.
}
	\item{cut.all.over.length}{If the chromosome is shorter than a position given in the interval data, an error is posted, except if this argument is \code{TRUE}, in which case the position is replaced by the length of the chromosome, or removed altogether. The default is \code{FALSE}.
}
  \item{ecdf.area.permut.number}{
The number of permutations performed to get the \emph{p-value} for the area between the ecdf for relative distance distribution and the straight line representing the uniform relative area for the random case. See below for explanation.
}
  \item{mean.distance.permut.number}{
The number of permutations to ascribe \emph{p-value} to minimal query-reference distance averaged over all query points.
}
  \item{jaccard.measure.permut.number}{
The number of permutations for Jaccard measure \emph{p-value} estimation. 
}
  \item{jaccard.permut.is.rearrangement}{
If \code{TRUE}, the permutations of the reference for the Jaccard test retain the lengths of all intervals and gaps in the query. All the permuted queries will mirror the original, so the \emph{p-value} is overestimated. If \code{FALSE} (the default), the permutation is a random resampling of starts of the query intervals.
}
  \item{awhole.space.name}{
The name of the pseudo-space that describes the overall genome statistics. Default is 'awhole'.
}
  \item{keep.distributions}{
It this is true, the procedure returns all points in th distributions calculated for comparison. This is useful for making figures. Default is \code{FALSE}.
}
  \item{representing.point.function}{
By default, the midpoint of each interval is used as the surrogate for the position of the interval. To force the program to use something other than the midpoint, define the function to use to return comparison points. The function must take the same parameters as the default \code{mitl} that returns the middle points. The function is to be passed as the \code{representing.point.function} parameter. The default for the parameter is:
\code{mitl<-function(start,end,chromosome.length,space){return ((as.integer(start)+as.integer(end))/2)}} 
}
  \item{query.representing.point.function}{
The same thing as the \code{representing.point.function}, but the representation point calculation is overloaded only for query intervals.
} 
  \item{reference.representing.point.function}{
The same thing as the \code{representing.point.function}, but the representation point calculation is overloaded only for query intervals.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The result is an instance of the \code{\link{GenometriCorrResult-class}} that describes the run
parameters in its \code{\link{GenometriCorrResult-class}} \code{@config} slot and that extends a
\code{\linkS4class{namedList}}
list (created originally with a \code{list()} call with the results of the run.
Each element of the list is also a list that describes results for a space (chromosome); one of
them is 'awhole' (or other \code{awhole.space.name} if given) that describes the genome awhole, all others are named the same as the chromosomes and describe the chromosomewise statistics.
The elements of the 'awhole' and chromosomewise lists are statistical measures and some datasets.
The statistical measures are described in the \code{\link{GenometriCorr}} package help.
For further explanation, see the MarkupsIndependence vignette. 

Below is the description of the values of the list returned for each chromosome. 
\item{query.population}{Query points used in the comparisons.}
\item{reference.population}{Reference points used in the comparisons.}
\item{relative.distances.ks.p.value}{\emph{p-value} for local independence obtained by the Kolmogorov-Smirnov test for relative distances. }
\item{relative.distances.ecdf.deviation.area.p.value}{\emph{p-value} for local independence obtained by the permutation test for relative distances. }
\item{relative.distances.ecdf.area.correlation}{Has the same sign with the relative distance-based local correlation. }
\item{projection.test.p.value}{\emph{p-value} for chromosome-scale independence obtained by the projection test.}
\item{projection.test.lower.tail}{If TRUE, projection test shows negative correlation, real overlap is less than the expectation.}
\item{scaled.absolute.min.distance.sum.p.value}{\emph{p-value} for chromosome-scale null hypothesis as obtained by the permutations of the query points and the mean of the distances to the two closest reference points.}
\item{scaled.absolute.min.distance.sum.lower.tail}{If TRUE, the query points are closer to the reference points than expected (the absolute distance is lower than the expectation). }
\item{query.reference.intersection}{Intersection of reference and query, in bases.}
\item{query.reference.union}{Union of reference and query, in bases.}
\item{jaccard.measure}{Jaccard measure of query and reference overlap.}
\item{jaccard.measure.p.value}{The permutation-based evaluation of the \emph{p-value} for the obtained Jaccard measure, given the null hypothesis of independence.}
\item{jaccard.measure.lower.tail}{If \code{TRUE}, then Jaccard measure is lower that the expectation (overlap less than expected)}


The additional values that are returned if \code{keep.distributions=TRUE}
\item{relative.distances.data}{The original relative distances}
\item{relative.distances.ecdf.deviation.area}{The real value of the ECDF deviation area to be compared with the permutation to obtain the p-value}
\item{relative.distances.ecdf.deviation.area.null.list}{The null distribution}
\item{projection.test}{List of three values: \code{reference.length} is length of a chromosome; \code{reference.coverage} is length of that chromosome covered by reference intervale, and \code{query.hits} is the number of query points that fall into the reference intervals.}
\item{absolute.min.distance.data}{The distribution of query-reference distances}
\item{absolute.inter.reference.distance.data}{The distribution of reference-reference distances}
\item{scaled.absolute.min.distance.sum}{The value of the sum (i.e. mean) of scaled absolute distances}
\item{scaled.absolute.min.distance.sum.null.list}{The null distribution for the scaled absolute distances}
\item{jaccard.measure.null.list}{The null distribution of Jaccard measures in permutations}

}
\references{
%% ~put references to the literature/web site here ~
\url{http://genometricorr.sourceforge.net/}
}

\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan}

%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	The \code{\link{GenometriCorr}} documentation and the \href{../doc/GenometriCorr.pdf}{GenometriCorr package vignette}.
}
\examples{

library('rtracklayer')
library('GenometriCorr')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

refseq<-transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
cpgis<-import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr"))
seqinfo(cpgis)<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqnames(seqinfo(cpgis))]

pn.area<-10
pn.dist<-10
pn.jacc<-10

cpgi_to_genes<-GenometriCorrelation(cpgis,refseq,chromosomes.to.proceed=c('chr1','chr2','chr3'),ecdf.area.permut.number=pn.area,mean.distance.permut.number=pn.dist,jaccard.measure.permut.number=pn.jacc,keep.distributions=FALSE,showProgressBar=FALSE)

print(cpgi_to_genes)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
