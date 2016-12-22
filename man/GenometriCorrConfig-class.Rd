% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2016 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
%               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
%               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
%
\name{GenometriCorrConfig-class}
\Rdversion{1.1}
\docType{class}
\alias{GenometriCorrConfig-class}
\alias{GenometriCorrConfig}
\alias{show,GenometriCorrConfig-method}

\title{Configuration class \code{"GenometriCorrConfig"}: a list of parameters for GenometriCorrelation, and input file names}

\description{An object of this class can be initialized by an \code{.ini} file; it is inherited form a \code{namedList}lthat contains the full structure of \code{.ini} file. It restores the structure by call to \code{print}. The \code{\link{GenometriCorr}} - specific call is \code{run.config(conf)} that runs \code{\link{GenometriCorrelation}} with parameters that are stored in the configuration file, as shown in \code{template-config.ini} file in the \code{extdata} directory of the \code{\link{GenometriCorr}} package. For a usage example, see \code{\link{GenometriCorrResult-class}}}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("GenometriCorrConfig", src)}.
The call creates the GenometriCorrConfig object from the file \code{src}.
The file is read as a list (\code{\linkS4class{namedList}}). Each section is a sub-list; each record is an element of its section list. For the string \code{a=b}, \code{a} is the name and \code{b} is the value; if the string does not contain \code{=}, it is read as \code{""}.

So, \preformatted{
[foo]
bar=2
tap
}

will be mapped to: \code{config$foo$bar} with value of \code{"2"} and \code{config$foo$tap} with value of \code{""}.
The values that are used to configure the \code{\link{GenometriCorrelation}} are described in the section about \code{run.config} method.
}

\section{Slots}{
  \describe{
    \item{\code{.Data}:}{Object of class \code{"list"} the list of values; the main content of the object.}
    \item{\code{src}:}{Object of class \code{"character"} the name of the file the object was initiated from}
    \item{\code{names}:}{Object of class \code{"character"} the names to address the mail list are names of the config file sections. It is inherited from "\linkS4class{namedList}" }
  }
}
\section{Extends}{
Class \code{"\linkS4class{namedList}"}, directly.
Class \code{"\linkS4class{list}"}, by class "namedList", distance 2.
Class \code{"\linkS4class{vector}"}, by class "namedList", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "GenometriCorrConfig")}: prints the \code{.ini} file text as a record or for use later as a \code{GenometriCorrConfig}.}
		}
}

\references{
%% ~put references to the literature/web site here ~
\url{http://genometricorr.sourceforge.net/}
}

\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan}

\note{
%%  ~~further notes~~
Here is an example of the \code{.ini} file content:
\preformatted{
[data]
query=UCSCcpgis_hg19.bed
query.format=bed
reference=UCSCrefseqgenes_hg19.bed
reference.format=bed
mapping=q.bed
do.mapping=FALSE

[chromosomes]
chr1
chr2
chr3

[chromosomes.length]
chr1=249250621
chr2=243199373
chr3=198022430
chr4=191154276
chr5=180915260
chr6=171115067
chr7=159138663
chrX=155270560
chr8=146364022
chr9=141213431
chr10=135534747
chr11=135006516
chr12=133851895
chr13=115169878
chr14=107349540
chr15=102531392
chr16=90354753
chr17=81195210
chr18=78077248
chr20=63025520
chrY=59373566
chr19=59128983
chr22=51304566
chr21=48129895
chrM=16571

[options]
add.chr.as.prefix=FALSE
awhole.only=FALSE
suppress.evaluated.length.warning=FALSE
cut.all.over.length=FALSE
keep.distributions=TRUE
showTkProgressBar=TRUE
showProgressBar=FALSE

[tests]
ecdf.area.permut.number=100
mean.distance.permut.number=100
jaccard.measure.permut.number=100
}

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
	The \code{\link{GenometriCorrResult-class}} class; the \code{\link{GenometriCorr}} documentation and the \href{../doc/GenometriCorr.pdf}{GenometriCorr package vignette}.
}

\keyword{classes}
