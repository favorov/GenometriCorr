% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2020 Alexander Favorov, Loris Mularoni, 
%               Yulia Medvedeva, Harris A. Jaffee, 
%               Ekaterina V. Zhuravleva, Veronica Busa,
%               Leslie M. Cope, Andrey A. Mironov, 
%               Vsevolod J. Makeev, Sarah J. Wheelan.
\name{GenometriCorr}
\alias{GenometriCorr}
\alias{GenometriCorr-package}
\docType{package}
\title{
GenometriCorr (Genometric Correlation) package 
}
\description{
\code{GenometriCorr} evaluates the spatial correlation between two types of features (experimental or annotated) in genomic coordinates. Several different approaches are implemented, each intended to evaluate a different biologically relevant spatial relationship. The genomic features are split into intervals and correlations are evaluated as described in the package vignette. 

Most standard genomic feature annotation file formats are accepted as input. Each feature type is given in a separate file, one used as a query and the other used as a reference. The statistical comparisons that are made are not symmetric; that is, they are sensitive to which file is the query and which is the reference, so we recommend running the comparisons twice, switching the status of the input files for the second run. This is intentional, as biological relationships can be asymmetric. Once the files are read, the annotations are stored as \code{\link{GRanges}} objects from \code{\link{GenomicRanges}} package).

There are two classes that are provided by the package. One of them, \code{\link{GenometriCorrConfig-class}} provide serialization into configuration file that describe a full run of the package that is intended to explore spatial correlations between two markups. The class also provide run.config method that open the input data files, read them, make all the calculations by applying the central \code{\link{GenometriCorrelation}} function to mapped or unmapped data (see below). The \code{\link{GenometriCorrelation}} function returns an instance of \code{\link{GenometriCorrResult-class}} that is based on a list with the results of the run and also contains the \code{\link{GenometriCorrConfig-class}} instance that describes the run as a slot. The \code{\link{GenometriCorrResult-class}} provides \code{\link{show}} method and two graphical representations: \code{\link{graphical.report}} and \code{\link{visualize}}. 

Below we describe the statistical comparisons implemented in this package. One set of features is termed the \dQuote{query} and the other, the \dQuote{reference}, throughout. 

\itemize{

	\item \code{scaled.absolute.min.distance.sum.p.value}
Query and reference intervals that are often separated by the same distance will have significant \emph{p-values}; the magnitude of the \emph{p-value} depends on the number of permutations performed. To determine whether a significant result indicates a small or large distance, we provide the \code{scaled.absolute.min.distance.sum.lower.tail}, a boolean that is \code{TRUE} when the actual absolute distances are smaller than expected and \code{FALSE} if they are higher.

	\item{relative.distances.ks.p.value}
This \emph{p-value} will be low when the query and reference intervals have similar spatial distributions. This test is hypersensitive and is most useful when the \emph{p-value} is high.

	\item{relative.distances.ecdf.deviation.area.p.value} 
If the query and reference reference intervals are closer or farther apart than expected, this \emph{p-value} will be lo
w.

	\item{relative.distances.ecdf.area.correlation}
This value will be negative when the query and reference intervals are anticorrelated and positive when they are positively correlated. This value has no relation to the \emph{p-value} of correlation.

	\item{projection.test.p.value}
If overlaps between query characteristic points (by default, midpoints) and the reference features occur less often than expected, the \code{projection.test.lower.tail} is \code{TRUE}; if they are more common than expected, it is \code{FALSE}. To measure the effect size, the observed to expected ratio is calculated as \code{projection.test.obs.to.exp}

	\item{jaccard.measure}
The Jaccard test compares the length of the union of all query and reference features with the length of the intersection of the query and reference features.
 
	\item{jaccard.measure.p.value}
This is the permutation-based evaluation of the \emph{p-value} for the Jaccard measure; the \code{jaccard.measure.lower.tail} is \code{TRUE} if there are fewer overlaps than expected (less overlap) ond \code{FALSE} otherwise.
}

When both query and reference features are restricted to genomic subsets, or when we want to compare their relationship only within smaller genomic intervals (e.g. genes), we can remap the intervals of interest onto pseudochromosomes and use the \code{\link{MapRangesToGenomicIntervals}} provided by the package. Each mapping result is a \code{\link{GRanges}} object. The pseudochromosomes can now be treated as usual by \code{\link{GenometriCorrelation}} to test the correlations of interest.


The package also provides the \code{\link{VisualiseTwoIRanges}} function that creates a very high-level graphic overview of a pair of annotations on a chromosome (space).

}
\details{
\tabular{ll}{
Package: \tab GenometriCorr\cr
Type: \tab Package\cr
Title: \tab Genometric Correlation package \cr
Version: \tab 1.1.23\cr
Date: \tab 2020-02-20\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
biocViews: \tab Annotation, Genetics, Infrastructure, DataRepresentation, Bioinformatics, StatisticalMethod\cr
URL: \tab http://genometricorr.sourceforge.net/ \cr
}}

\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Veronica Busa, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan}

\references{
\url{http://genometricorr.sourceforge.net/}
}
\keyword{ package }
\seealso{ See the GenometriCorr package vignette.}

\examples{
library('rtracklayer')
library('GenometriCorr')

library('rtracklayer')
library('TxDb.Hsapiens.UCSC.hg19.knownGene')

refseq<-transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)

cpgis<-import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr"))
seqinfo(cpgis)<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqnames(seqinfo(cpgis))]

permut.number<-0

#permut.number=0 means all the permutations are off
#permut.number is the common default for ecdf.area.permut.number, mean.distance.permut.number, and jaccard.measure.permut.number
#these three can be set separetely; explicit set overrloads the default


cpgi_to_genes<-GenometriCorrelation(cpgis,refseq,chromosomes.to.proceed=c('chr1','chr2','chr3'),permut.number=permut.number,keep.distributions=FALSE,showProgressBar=FALSE)

print(cpgi_to_genes)

VisualiseTwoIRanges(
	ranges(cpgis[seqnames(cpgis)=='chr1']),
	ranges(refseq[seqnames(refseq)=='chr1']),
	nameA='CpG Islands',nameB='RefSeq Genes',
	chrom_length=seqlengths(TxDb.Hsapiens.UCSC.hg19.knownGene)['chr1'],
	title="CpGIslands and RefGenes on chr1 of Hg19 animal")

#mapping example, same as in the vignette
population<-1000

chromo.length<-c(3000000)

names(chromo.length)<-c('the_chromosome')

rquery<-GRanges(ranges=IRanges(start=runif(population,1000001,2000000-9),width=c(10)),seqnames='the_chromosome')

rref<-GRanges(ranges=IRanges(start=runif(population,1000001,2000000-9),width=c(10)),seqnames='the_chromosome')

#create two features, they are randomly scattered in 1 000 000...2 000 000

unmapped_result<-GenometriCorrelation(rquery,rref,chromosomes.length=chromo.length,permut.number=permut.number,keep.distributions=FALSE,showProgressBar=FALSE)

#correlate them on the whole chromosome: 1...3 000 000

cat('Unmapped result:\n')
print(unmapped_result)

map_space<-GRanges(ranges=IRanges(start=c(1000001),end=c(2000000)),seqname='the_chromosome')

mapped_rquery<-MapRangesToGenomicIntervals(what.to.map=rquery,where.to.map=map_space)

mapped_rref<-MapRangesToGenomicIntervals(what.to.map=rref,where.to.map=map_space)

#map them into 1 000 001...2 000 000

mapped_result<-GenometriCorrelation(mapped_rquery,mapped_rref,permut.number=permut.number,keep.distributions=FALSE,showProgressBar=FALSE)

#then, correlate again

cat('Mapped result:\n')
print(mapped_result)

}
