% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2018 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
%               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
%               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
%
\name{readTableToIRanges}
\alias{readTableToIRanges}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
readTableToIRanges
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
readTableToIRanges.R reads a variety of GFF- or BED-like (text table) formats to an IRanges or a RangedData object.
It tries to interpret the header line if \code{header} parameter is TRUE.
}
\usage{
readTableToIRanges(file = NA, space = NA,
    start = NA, end = NA, width = NA, 
    sep = "", skip = 0, comment.char = "#", 
    header = F, autofeatures = F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
		\code{file} is the name of file or URL; it is passed to read.table so it was be anything that read.table allows as a file parameter.    
	}
  \item{header}{
		\code{header} indicates whether the first non-comment row should be interpreted as a header row. If header is TRUE, each column is named by the its label in the header raw. Also, if a header is given, the function tries to fit the given column names to the expected chrom, end, start, etc. The default is \code{FALSE}, if there is a header is the data, the function should be called with \code{header=TRUE} or with \code{skip=an integer} parameter.  
	}
  \item{space}{
		\code{space} is the number (1-based) or the name in the header of the column containing the chromosome name in the RangedData notation (the 'space'). If space is NA or is not given and not found in the header, the data are read to IRanges object with exceptions as listed below.
	}
  \item{start}{
		\code{start} gives the column number or column name for the interval start points
	}
  \item{end}{
		\code{end} gives the column number or column name for the interval end points
	}
  \item{width}{
		\code{width} gives the column number or column name for the interval start points; both \code{end} and \code{width} cannot be given
	}
  \item{sep}{
		\code{sep} is the field separator used in the file; "" means a whitespace-delimited table
	}
  \item{skip}{
		\code{skip} is the number of lines to be ignored at the beginning of the file
	}
  \item{comment.char}{
		\code{comment.char} gives the prefix of any lines that should be ignored
	}
  \item{autofeatures}{
		If \code{header} is true, then \code{autofeatures} keeps track of all columns not used by the function and adds them to the result 
	}
  \item{\dots}{
		\code{\dots} contains parameters in this form: feature_name=col_number or feature_name=column name, where the name or column is not one specified by other parameters. Each of these pairs will add the \code{feature_name} to the data in the specified column when it is returned with the result.
	}
}

%%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
Returns an IRanges or RangedData object. The type depends on parameters and the data file. Data that span chromosomes will have \code{space}, or \code{autofeatures=TRUE}, or \code{\dots}, to break the given intervals into subsets based on \code{spaces}. If \code{space} is not given or is NA, then all intervals are taken to be on the same chromosome (in the same space), and the name of the space is the file name.
}
\references{
%% ~put references to the literature/web site here ~
\url{http://genometricorr.sourceforge.net/}
}

\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan}

\note{
%%  ~~further notes~~
If the data you have is standard UCSC-like .bed, do not use this fuction, use \code{import} from \code{rtracklayer}. \code{import} correctly accounts for the difference in base numbering: UCSC assumes 0-based interval starts and 1-based interval ends while IRanges assumes that all data are 1-based.
The \code{readTableToIRanges} function is necessary when you work with nonstandard data in table form.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
	The \code{\link{GenometriCorr}} documentation and the \href{../doc/GenometriCorr.pdf}{GenometriCorr package vignette}.
	\code{\link[rtracklayer:export]{Import}} import function from \code{rtracklayer} 
}

\examples{
library('rtracklayer')
library("GenometriCorr")

cpgi_data<-system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr")

cpgis_imported<-as(import(cpgi_data),'RangedData');

cpgis_read<-GenometriCorr:::readTableToIRanges(cpgi_data,comment.char='$',header=TRUE);

#the starts the BED file are 0-based,
#import() raised them for cpgis_imported
#here we do it manually for cpgis_read
start(cpgis_read)<-start(cpgis_read)+1;

#rtraclayer added flowery metadata
metadata(ranges(cpgis_imported))<-list();

if (identical(cpgis_read,cpgis_imported))
{
	cat("The two results are identical.\n")
}
#now, they are identical


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
\keyword{ files }% __ONLY ONE__ keyword per line
