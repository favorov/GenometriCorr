% GenometriCorr (Genometric Correlation): an R package for spatial correlation of genome-wide interval datasets
% (c) 2010-2011 Alexander Favorov, Loris Mularoni, Leslie M. Cope, Andrey A. Mironov  
%              Yulia Medvedeva, Vsevolod J. Makeev, Sarah J. Wheelan
%
% favorov@sensi.org
%$Id: MapRangesToGenomicIntervals.Rd 1715 2012-04-18 15:06:29Z favorov $


\name{MapRangesToGenomicIntervals}
\alias{MapRangesToGenomicIntervals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
MapRangesToGenomicIntervals
}
\description{
MapRangesToGenomicIntervals is a function that maps a set of ranges represented by an \code{RangedData} or by \code{GRanges} object into a set of intervals (subgenome) (e.g. RefSeq) that is also represented by an \code{RandedData} or \code{GRanges}.

The result is a \code{GRanges}. Each space in the output value corresponds to an interval from \code{where.to.map}, and the items in \code{what.to.map} have been transferred to the new coordinate system in each interval.

The function is extremely useful if two features are spatially restricted in the genome. See more in the \href{../doc/GenometriCorr.pdf}{GenometriCorr package vignette}.  

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MapRangesToGenomicIntervals(
	where.to.map, what.to.map, 
	chromosomes.to.proceed=NA, 
	unmapped.chromosome.warning=TRUE, 
	unmapped.range.warning=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{where.to.map}{
		The set of genomic intervals we map to.
	}
  \item{what.to.map}{
		The set of ranges that we map.
	}
  \item{chromosomes.to.proceed}{
		The default set of chromosomes to map is the intersection of the chromosomes in where.to.map and what.to.map. If we want to restrict the set, we can do it with this parameter.
	}
	\item{unmapped.chromosome.warning}{
		For each chromosome that is represented in \code{what.to.map} and that is included in \code{chromosomes.to.proceed} if it is given and that is not represented in \code{where.to.map}, a warning is generated if \code{unmapped.chromosome.warning} is \code{TRUE}. The default is \code{TRUE}.
	}
	\item{unmapped.range.warning}{
		If a genomic region in \code{what.to.map} cannot be mapped and is contained in \code{chromosomes.to.proceed}, a warning is generated if \code{unmapped.range.warning} is \code{TRUE}. The default is \code{FALSE}.  
	}
}
\details{
We assume here that all of the ranges are 1-based and end-inclusive (end is 1-based, too).
}

\value{
The result is a \code{GRanges} object.  Each space in the output is an interval from \code{where.to.map}. The ranges of the value are those from \code{what.to.map} with their coordinates mapped to the new spaces.

The value contains all those ranges from \code{what.to.map} that fall into an interval from \code{where.to.map} with their coords relative to \code{where.to.map} the interval start; the space (chromosome) names are now like \code{chr1_1000} that means that the range is mapped on the \code{where.to.map} interval that strats at position 1000 of chromosome chr1.
}

\references{
%% ~put references to the literature/web site here ~
\url{http://genometricorr.sourceforge.net/}
}
\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Leslie M. Cope, Andrey A. Mironov,  Yulia Medvedeva, Vsevolod J. Makeev, Sarah J. Wheelan}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	The \code{\link{GenometriCorr}} documentation and the \href{../doc/GenometriCorr.pdf}{GenometriCorr package vingnette}.
}
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\examples{
library('GenometriCorr')

intervals<-RangedData(ranges=IRanges(c(1,10001,1,10001),width=c(1000)),space=c('chr1','chr1','chr2','chr2'))
ranges=RangedData(ranges=IRanges(c(10,110,10101,100000,500,550,1055),width=c(10)),space=c(rep('chr1',4),rep('chr2',3)))

mapped<-MapRangesToGenomicIntervals(where.to.map=intervals,what.to.map=ranges)

#the result is:
#GRanges with 5 ranges and 0 elementMetadata values
#      seqnames     ranges strand |
#         <Rle>  <IRanges>  <Rle> |
#[1]     chr1_1 [ 10,  19]      * |
#[2]     chr1_1 [110, 119]      * |
#[3] chr1_10001 [101, 110]      * |
#[4]     chr2_1 [500, 509]      * |
#[5]     chr2_1 [550, 559]      * |
#
#seqlengths
#     chr1_1 chr1_10001     chr2_1 chr2_10001
#       1000       1000       1000       1000

}
\keyword{manip}% __ONLY ONE__ keyword per line
