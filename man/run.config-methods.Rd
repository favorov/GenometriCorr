% GenometriCorrelation project evaluating two interval markups genomewide independence. 
% (c) 2010-2020 Alexander Favorov, Loris Mularoni, 
%               Yulia Medvedeva, Harris A. Jaffee, 
%               Ekaterina V. Zhuravleva, Veronica Busa,
%               Leslie M. Cope, Andrey A. Mironov, 
%               Vsevolod J. Makeev, Sarah J. Wheelan.
%
\name{run.config-methods}
\docType{methods}
\alias{run.config}
\alias{run.config-methods}
\alias{run.config,GenometriCorrConfig-method}
\title{Method for Function \code{run.config} in Package \pkg{GenometriCorr} }
\description{
The \code{run.config} function runs on a \code{"\linkS4class{GenometriCorrConfig}"} class object; it runs the \code{\link{GenometriCorrelation}} function with the parameters described in the configuration object. The resulting object of \code{"\linkS4class{GenometriCorrResult}"} has the original options in its \code{config} slot.
}
\section{Methods}{
\describe{
	\item{\code{signature(conf = "GenometriCorrConfig")}}{
		\code{run.config(conf,query=NA,reference=NA,query.format=NA,reference.format=NA,mapping=NA,mapping.format=NA,do.mapping=NA)}\cr
		runs \code{\link{GenometriCorrelation}}.\cr
		\code{conf} is the \code{"\linkS4class{GenometriCorrConfig}"} class. Its slots correspond to the configuration file sections. Below is the list of the sections and the variables of the file and the corresponding parameters of the \link{GenometriCorrelation}} call. 
		\tabular{rll}{
		Section \tab Value \tab Role\cr
		\code{[chromosomes]} \tab \tab this section contains the list of chromosomes, it is passed as \code{chromosomes.to.proceed}\cr
		\code{[chromosomes.length]}\tab \code{name=length} \tab this length will be passed in the \code{chromosomes.length} parameter in the \code{['name']} element\cr
		\code{[options]} \tab \code{add.chr.as.prefix} \tab \code{add.chr.as.prefix}\cr
		\code{[options]} \tab \code{awhole.only} \tab \code{awhole.only}\cr
		\code{[options]} \tab \code{suppress.evaluated.length.warning} \tab \code{suppress.evaluated.length.warning}\cr
		\code{[options]} \tab \code{keep.distributions} \tab \code{keep.distributions}\cr
		\code{[options]} \tab \code{showTkProgressBar} \tab \code{showTkProgressBar}\cr
		\code{[options]} \tab \code{showProgressBar} \tab \code{showProgressBar}\cr
		\code{[tests]} \tab \code{permut.number} \tab \code{permut.number}\cr
		\code{[tests]} \tab \code{ecdf.area.permut.number} \tab \code{ecdf.area.permut.number}\cr
		\code{[tests]} \tab \code{mean.distance.permut.number} \tab \code{mean.distance.permut.number}\cr
		\code{[tests]} \tab \code{jaccard.measure.permut.number} \tab \code{jaccard.measure.permut.number}\cr
		\code{[tests]} \tab \code{random.seed} \tab value to seed the random generator (for reproducibility)
		}
		Below are values that describe the input data; they can be overriden by the \code{run.config} call parameters.
		\tabular{rll}{
		Section \tab Value \tab Role\cr
		\code{[data]} \tab \code{query} \tab name of the query file\cr
		\code{[data]} \tab \code{query.format} \tab format of the file or \code{'R.variable.name'} \cr
		\code{[data]} \tab \code{reference} \tab name of the reference file\cr
		\code{[data]} \tab \code{reference.format} \tab format of the file or \code{'R.variable.name'} \cr
		\code{[data]} \tab \code{mapping} \tab name of the mapping file\cr
		\code{[data]} \tab \code{mapping.format} \tab format of the mapping file or \code{'R.variable.name'} \cr
		\code{[data]} \tab \code{do.mapping} \tab logical; whether to do mapping
		}
		Each of the data (\code{query}, \code{reference}, \code{mapping}) parameters may instead be data of one of the accepted formats (\code{\link{GRanges}}). In this case, the corresponding format value is omitted and nothing is written to the \code{config} slot of the result, and the data are used for the function calls. If the format parameter is \code{'R.variable.name'}, the input is taken from that variable. The calculation is the same as for \code{\link{GRanges}} parameter, but the name of the variable and the \code{'R.variable.name'} are saved in \code{config} slot of the result.

		Other formats are interpreted as \code{format} parameters for the {\code{\link[rtracklayer:export]{import}}} function from \code{rtracklayer} package. 
}}

\author{Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Veronica Busa,Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan}

\note{See an example at \code{\linkS4class{GenometriCorrResult}} help.}

\keyword{methods}
