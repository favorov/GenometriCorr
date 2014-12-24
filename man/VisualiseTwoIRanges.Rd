% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2014 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
%               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
%               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
%
\name{VisualiseTwoIRanges}
\alias{VisualiseTwoIRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
VisualiseTwoIRanges
}
\description{
VisualiseTwoIRanges is a function that displays the intervals within two IRanges, one above the other, in different colors, along a chromosome or subset of a chromosome. The intent is to show large-scale relationships between the two IRanges.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VisualiseTwoIRanges(irA, irB,
    start = 1, end = NA,
    nameA = "RangesA", nameB = "RangesB", 
    chrom_length = NA, title = NA,
    pdf = NULL, close.device=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{irA}{

	}
  \item{irB}{
		The two IRanges objects to visualise
	}
  \item{start}{

	}
  \item{end}{
		Boundaries of the area to be pictured
	}
  \item{nameA}{

	}
  \item{nameB}{
		The names of the IRanges data to use as the plot legend	and for the file name if \code{pdf=""} is given as parameter.}
  \item{chrom_length}{
		The length of the chromosome spanned by the IRanges
	}
	\item{title}{Title, printed at the top of the plot.
	}

   \item{pdf}{
		Name of a file to which the image should be written.
		If \code{pdf=""} the filename is constructed from \code{nameA} and \code{nameB}.
		If \code{NULL}, no pdf is opened and an x11 window is raised up as by \code{plot}.
		The suffix ".pdf" will be appended if not included.
	}
   \item{close.device}{
		Whether to close the plot device after writing the image.
		A \code{FALSE} setting allows multiple images to be written
		to the same pdf file, but then eventually closing the
		device is up to the user.
		Default is \code{NULL} that means: close an x11, do npt close a pdf.
	}
}
%\details{
%}
%\value{
%}
\references{
%% ~put references to the literature/web site here ~
\url{http://genometricorr.sourceforge.net/}
}
\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Leslie M. Cope, Andrey A. Mironov,  Yulia Medvedeva, Vsevolod J. Makeev, Sarah J. Wheelan}
%\note{
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
	The \code{\link{GenometriCorr}} documentation and the {GenometriCorr package vignette}.
}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\examples{

library('rtracklayer')
library('GenometriCorr')

cpgis<-as(import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr")),'RangedData');

refseq<-as(import(system.file("extdata", "UCSCrefseqgenes_hg19.bed", package = "GenometriCorr")),'RangedData');


human.chrom.length<-c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,59373566,155270560)


names(human.chrom.length)<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX')

VisualiseTwoIRanges(cpgis['chr1']$ranges, refseq['chr1']$ranges,
	nameA='CpG Islands', nameB='RefSeq Genes',
	chrom_length=human.chrom.length[['chr1']],
	title="CpGIslands and RefGenes on chr1, hg19",
	pdf='CpGi_vs_RefSeq_genes_chr1_hg19')

}
\keyword{ hplot }% __ONLY ONE__ keyword per line
