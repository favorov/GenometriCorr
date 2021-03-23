# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2020 Alexander Favorov, Loris Mularoni, 
#               Yulia Medvedeva, Harris A. Jaffee, 
#               Ekaterina V. Zhuravleva, Veronica Busa,
#               Leslie M. Cope, Andrey A. Mironov, 
#               Vsevolod J. Makeev, Sarah J. Wheelan.
#
# GenomertiCorrelation is the main function of the package

#' @importFrom gtools mixedsort
#' @importFrom stringr str_detect str_c
#' @importFrom tcltk tkProgressBar setTkProgressBar getTkProgressBar
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqlengths
#' @importFrom stats ecdf integrate ks.test pbinom punif runif 
#' @importFrom utils getTxtProgressBar head packageDescription read.table setTxtProgressBar tail txtProgressBar 
#' @import BiocGenerics magrittr GenomicRanges GenomicFeatures methods plyranges

epsilon=1e-6
integr_rel_tol=0.01

GenometricCorrelation <- function(...)
{
   .Deprecated("GenometriCorrelation")
   GenometriCorrelation(...)
}

#we are going to use it to convert chromosome spaces to strings/string lists
listtoString<-function(mylist)
{
	if (! is.vector(mylist) && !is.factor(mylist) ) return (c())
	sapply(mylist,function(x){toString(x)})
}

mitl<-function(start,end,chromosome.length,space){return ((as.integer(start)+as.integer(end))/2)}
#middle


add.chr.prefix.to.names<-function(namelist)
#adds 'chr' prefix to all names if they are not 'chr'-prefixed already;
#lowercases the 'CHR','ChR', etc prefixes if any instead of adding 'chr'
{
	if (inherits(namelist,"GRanges"))
	{	
		sinf<-namelist@seqinfo
		sinf@seqnames<-add.chr.prefix.to.names(sinf@seqnames)
		namelist<-GRanges(
			ranges=namelist@ranges,
			seqnames=add.chr.prefix.to.names(namelist@seqnames),
			seqinfo=sinf
		)
		return(namelist)
	}
	if (is.factor(namelist))
	{
		levels(namelist)<-add.chr.prefix.to.names(levels(namelist))
		return(namelist)
	}
	
	result<-list()

	for ( name in namelist )
	{
		if (is.vector(name) && length(name)>1)
			result[[length(result)+1]]<-add.chr.prefix.to.names(name)
			#tree recursion
		else
		{
			if (!is.character(name)) stop('Non-character element passed to add.chr.prefix.to.names.\n')
			if (tolower(substr(name,1,3))	=='chr')
				result<-c(result,paste('chr',substring(name,4),sep=''))
			else
				result<-c(result,paste('chr',name,sep=''))
		}
	}
	if (!is.list(namelist))
		return(as.character(result))
	else 
		return(result)
}


#' GenometriCorrelation
#'
#' The GenometriCorrelation function compares two interval annotations on a chromosome, set of chromosomes or on an entire genome, and performs various statistical tests to determine whether the two interval sets are independent or are positioned nonrandomly with respect to each other. For a complete description of the tests, refer to the the package vignette.
#'
#' @aliases GenometricCorrelation
#' @param query GRanges object that contains the query interval set coordinates and spaces (generally, chromosomes). The \code{GenometriCorrelation} function tests whether it is positioned independently relative to the \code{reference} interval set.
#' @param reference The \code{GRanges} object that contains the reference interval set coordinates and spaces. The \code{GenometriCorrelation} function tests whether \code{query} is positioned independently relative to \code{reference}.
#' @param chromosomes.to.proceed This vector of strings contains the names of spaces (chromosomes) to analyze. If both \code{query} and \code{reference} are \code{GRanges} if the parameter is not given (its default is c()), the initial list of spaces to proceed is the intersection of the space lists of query and reference names, and a warning is rasied if the lists are different. if \code{chromosomes.to.proceed} is given, it restricts the list so that only those names from the intersection that are in \code{chromosomes.to.proceed} are analyzed. If it is a subset of the \code{query} and \code{reference} chromosone lists intersection, the warning that the lists are different is not raised.
#' @param chromosomes.to.include.in.awhole ia a vector of strings contains the names of spaces to be included in the overall (awhole) statistics. Its default is c(), meaning that all the analysed genes are included.
#' @param chromosomes.to.exclude.from.awhole is a list of chromosomes (spaces) to be excluded from the overall statistics.
#' @param add.chr.as.prefix deals with the chr chromosome name prefix. The correlation is only performed on chromosomes that have exactly the same name, so by default, a chromosome named chr1 will not be considered the same chromosome as one simply labeled 1. This argument is provided so that if the chromosome names in the \code{query}, \code{reference} and the \code{chromosomes*} parameters, or the names of the chromosomes in \code{chromosome.length}) have no prefix \code{'chr'} (lower-, upper-, or mixed case), this prefix will be added, and all strings matching \code{'chr'} in uppercase or partly uppercase are changed to lowercase. This way, comparisons can be made even if the chromosome names differ by any variation of the word \code{'chr'}. The default is FALSE.
#' @param awhole.only If \code{FALSE}, all the considered chromosomes statistics and the summary (awhole) information are returned (default). If \code{TRUE}, only the summary (awhole) information is returned. Useful for pipelines. The \code{TRUE} value is default that is posted by \code{run.config} if the data is processed through mapping.
#' @param map.to.half Some of the tests we use are besed on distances between a query point and the closest reference point. If map.to.half is TRUE (default) we look for the closest reference point upstream or downstream, if it is \code{FALSE}, we look only downstream (left). Useful if you are mapping to intervals and want to preserve directionality.
#' @param showProgressBar Toggle the text progress bar. Default is \code{TRUE}.
#' @param showTkProgressBar Toggle the Tk progress bar. If it is \code{TRUE} but the Tcl/Tk is not loadable (e.g. \code{require('tcltk')} returns \code{FALSE}), the parameter is turned to \code{FALSE}. Default is \code{FALSE}.
#' @param chromosomes.length A vector of lengths of the chromosomes to be tested; each chromosome is given a name and a numerical length. The order of the chromosomes does not matter.
#' @param suppress.evaluated.length.warning If there is a chromosome that is included in the evaluation but its length is not given in \code{chromosome.length}, the length of the chromosome is calculated as the maximal position mentioned in the data. If this parameter is \code{FALSE}, this will generate a warning, which is suppressed when the parameter is \code{TRUE}. The default is \code{FALSE}.
#' @param permut.number is the common default for \code{ecdf.area.permut.number}, \code{mean.distance.permut.number}, and \code{jaccard.measure.permut.number}. \code{permut.number=0} defaults all the permutations to off.
#' @param ecdf.area.permut.number The number of permutations performed to get the \emph{p-value} for the area between the ecdf for relative distance distribution and the straight line representing the uniform relative area for the random case.
#' @param mean.distance.permut.number The number of permutations to ascribe \emph{p-value} to minimal query-reference distance averaged over all query points.
#' @param jaccard.measure.permut.number The number of permutations for Jaccard measure \emph{p-value} estimation. 
#' @param jaccard.permut.is.rearrangement If \code{TRUE}, the permutations of the reference for the Jaccard test retain the lengths of all intervals and gaps in the query. All the permuted queries will mirror the original, so the \emph{p-value} is overestimated. If \code{FALSE} (the default), the permutation is a random resampling of starts of the query intervals.
#' @param alternative a character string specifying the alternative hypothesis, must be one of \code{"two.sided"} (default), \code{"attraction"} or \code{"repulsion"}. You can specify just the initial letter.
#' @param awhole.space.name The name of the pseudo-space that describes the overall genome statistics. Default is 'awhole'.
#' @param keep.distributions It this is true, the procedure returns all points in th distributions calculated for comparison. This is useful for making figures. Default is \code{FALSE}.
#' @param representing.point.function By default, the midpoint of each interval is used as the surrogate for the position of the interval. To force the program to use something other than the midpoint, define the function to use to return comparison points. The function must take the same parameters as the default \code{mitl} that returns the middle points. The function is to be passed as the \code{representing.point.function} parameter. The default for the parameter is: \code{mitl<-function(start,end,chromosome.length,space){return ((as.integer(start)+as.integer(end))/2)}} 
#' @param query.representing.point.function The same thing as the \code{representing.point.function}, but the representation point calculation is overloaded only for query intervals.
#' @param reference.representing.point.function The same thing as the \code{representing.point.function}, but the representation point calculation is overloaded only for query intervals.
#' @param supress.evaluated.length.warning It was a typo for \code{suppress.evaluated.length.warning}. Now obsoleted, use \code{suppress.evaluated.length.warning}
#' @return The result is an instance of the \code{\link{GenometriCorrResult-class}} that describes the run parameters in its \code{\link{GenometriCorrResult-class}} \code{@config} slot and that extends a \code{\linkS4class{namedList}} list (created originally with a \code{list()} call with the results of the run.  Each element of the list is also a list that describes results for a space (chromosome); one of them is 'awhole' (or other \code{awhole.space.name} if given) that describes the genome awhole, all others are named the same as the chromosomes and describe the chromosomewise statistics.  The elements of the 'awhole' and chromosomewise lists are statistical measures and some datasets. The statistical measures are described in the \code{\link{GenometriCorr}} package help.  For further explanation, see the the package vignette. 
#'
#' Below is the description of the values of the list returned for each chromosome. 
#' \item{query.population}{Query points used in the comparisons.}
#' \item{reference.population}{Reference points used in the comparisons.}
#' \item{alternative}{Shows the value of the \code{alternative} parameter passed.} 
#' \item{relative.distances.ks.p.value}{\emph{p-value} for local independence obtained by the Kolmogorov-Smirnov test for relative distances. }
#' \item{relative.distances.ecdf.deviation.area.p.value}{\emph{p-value} for local independence obtained by the permutation test for relative distances. }
#' \item{relative.distances.ecdf.deviation.area.test.direction}{"attraction" or "repulsion". If the \code{alternative} parameter is "two.sided" (default), if shows the direction of the effect; if the parameter is {"attraction" or "repulsion"} it is not shown}
#' \item{relative.distances.ecdf.area.correlation}{Has the same sign with the relative distance-based local correlation. }
#' \item{projection.test.p.value}{\emph{p-value} for chromosome-scale independence obtained by the projection test. }
#' \item{projection.test.direction}{"attraction" or "repulsion". If the \code{alternative} parameter is "two.sided" (default), if shows the direction of the effect; if the parameter is {"attraction" or "repulsion"} it is not shown}
#' \item{projection.test.obs.to.exp}{To measure the effect size, the observed to expected ratio for the projection test statistics that is the number of query characteristic points (by default, midpoints) that fell into a reference features.}
#' \item{scaled.absolute.min.distance.sum.p.value}{\emph{p-value} for chromosome-scale null hypothesis as obtained by the permutations of the query points and the mean of the distances to the two closest reference points.}
#' \item{scaled.absolute.min.distance.sum.test.direction}{"attraction" or "repulsion". If the \code{alternative} parameter is "two.sided" (default), if shows the direction of the effect; if the parameter is {"attraction" or "repulsion"} it is not shown}
#' \item{query.reference.intersection}{Intersection of reference and query, in bases.}
#' \item{query.reference.union}{Union of reference and query, in bases.}
#' \item{jaccard.measure}{Jaccard measure of query and reference overlap.}
#' \item{jaccard.measure.p.value}{The permutation-based evaluation of the \emph{p-value} for the obtained Jaccard measure, given the null hypothesis of independence.}
#' \item{jaccard.measure.test.direction}{"attraction" or "repulsion". If the \code{alternative} parameter is "two.sided" (default), if shows the direction of the effect; if the parameter is {"attraction" or "repulsion"} it is not shown}
#'
#' The additional values that are returned if \code{keep.distributions=TRUE}
#' \item{relative.distances.data}{The original relative distances}
#' \item{relative.distances.ecdf.deviation.area}{The real value of the ECDF deviation area to be compared with the permutation to obtain the p-value}
#' \item{relative.distances.ecdf.deviation.area.null.list}{The null distribution}
#' \item{projection.test}{List of three values: \code{space.length} is length of a chromosome; \code{reference.coverage} is length of that chromosome covered by reference intervale, and \code{query.hits} is the number of query points that fall into the reference intervals.}
#' \item{absolute.min.distance.data}{The distribution of query-reference distances}
#' \item{absolute.inter.reference.distance.data}{The distribution of reference-reference distances}
#' \item{scaled.absolute.min.distance.sum}{The value of the sum (i.e. mean) of scaled absolute distances}
#' \item{scaled.absolute.min.distance.sum.null.list}{The null distribution for the scaled absolute distances}
#' \item{jaccard.measure.null.list}{The null distribution of Jaccard measures in permutations}
#' @author Alexander Favorov \email{favorov@@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Veronica Busa, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#' @references \href{http://genometricorr.sourceforge.net/}{GenometriCorr home}
#' @seealso The \code{\link{GenometriCorr}} documentation and vignette.
#' @examples
#' library('rtracklayer')
#' library('GenometriCorr')
#' library('TxDb.Hsapiens.UCSC.hg19.knownGene')
#' 
#' refseq<-transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' cpgis<-import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr"))
#' seqinfo(cpgis)<-seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[seqnames(seqinfo(cpgis))]
#' 
#' pn.area<-10
#' pn.dist<-10
#' pn.jacc<-10
#' 
#' cpgi_to_genes<-GenometriCorrelation(cpgis,refseq,chromosomes.to.proceed=c('chr1','chr2','chr3'),ecdf.area.permut.number=pn.area,mean.distance.permut.number=pn.dist,jaccard.measure.permut.number=pn.jacc,keep.distributions=FALSE,showProgressBar=FALSE)
#' 
#' print(cpgi_to_genes)

#' @keywords multivariate
#' @export
GenometriCorrelation <- function(
	query,reference,
	chromosomes.to.proceed=c(),
	chromosomes.to.include.in.awhole=c(),
	chromosomes.to.exclude.from.awhole=c(),
	add.chr.as.prefix=FALSE,
	awhole.only=FALSE,
	map.to.half=TRUE,
	showProgressBar=TRUE,
	showTkProgressBar=FALSE,
	chromosomes.length=c(),
	suppress.evaluated.length.warning=FALSE,
	permut.number=100,
	ecdf.area.permut.number=permut.number,
	mean.distance.permut.number=permut.number,
	jaccard.measure.permut.number=permut.number,
	jaccard.permut.is.rearrangement=FALSE,
	alternative='two.sided',
	awhole.space.name="awhole",
	keep.distributions=FALSE,
	representing.point.function=mitl,
	query.representing.point.function=representing.point.function,
	reference.representing.point.function=representing.point.function,
	supress.evaluated.length.warning
	)
{

	#fixing typo
	if('supress.evaluated.length.warning' %in% names(as.list(match.call())))
	{
		warning("Please use \'suppress.evaluated.length.warning\' instead of \'supress.evaluated.length.warning\'")
		if(supress.evaluated.length.warning==TRUE) 
		{
			suppress.evaluated.length.warning=TRUE
		}
	}
	#
	nameQ<-"query"
	nameR<-"reference"	# HJ fixed typo
	if (!is.object(query))  
		stop("The thing given as first (",nameQ,") range argument to\n  GenomertiCorrelation is not an object!")
	if (!is.object(reference))  
		stop("The thing given as second (",nameR,") range argument to\n  GenomertiCorrelation is not an object!")
	#they are both objects if we are here


	if (!(inherits(query,"GenomicRanges")))
		stop("The thing given as first (",nameQ,") range argument to GenomertiCorrelation is\n  not an GRanges object!")
	if (!(inherits(reference,"GenomicRanges")))
		stop("The thing given as second (",nameR,") range argument to GenomertiCorrelation is\n  not an GRanges object!")

	#they both are  GRanges if we are here

	if (!setequal(seqlevels(query),seqlevels(reference))) {
		common_seqinfo<-intersect(seqinfo(query),seqinfo(reference))
		common_seqs<-seqnames(common_seqinfo)
		if ( length(chromosomes.to.proceed)==0 || !all(chromosomes.to.proceed %in% common_seqs)) {
			warning("Query and referance has different chromosome lists.")
		}
		query<-plyranges::filter(query,seqnames %in% common_seqs)
		reference<-plyranges::filter(reference, seqnames %in% common_seqs)
		#actually, we just set seqinfo to common_seqinfo, but....
		query<-GRanges(seqnames=as.character(query@seqnames),ranges=query@ranges,
					strand=query@strand,mcols=mcols(query),seqinfo = common_seqinfo)
		reference<-GRanges(seqnames=as.character(reference@seqnames),ranges=reference@ranges,
					strand=reference@strand,mcols=mcols(reference),seqinfo = common_seqinfo)
	} else {
		common_seqs<-seqlevels(query) #they are equal, so quary and reference is the same
	}

	if ( length(chromosomes.to.proceed)>0 ){chromosomes.to.proceed<-intersect(common_seqs,chromosomes.to.proceed)}
	else {
		chromosomes.to.proceed<-common_seqs
	}
	#now, we see only sequences that are in both annotations and in chrosomes.to.proceed


	if (add.chr.as.prefix)
	{
		chromosomes.to.proceed<-add.chr.prefix.to.names(chromosomes.to.proceed)
		chromosomes.to.exclude.from.awhole<-add.chr.prefix.to.names(chromosomes.to.exclude.from.awhole)
		chromosomes.to.include.in.awhole<-add.chr.prefix.to.names(chromosomes.to.include.in.awhole)
		names(chromosomes.length)<-add.chr.prefix.to.names(names(chromosomes.length))
		query<-add.chr.prefix.to.names(query)
		reference<-add.chr.prefix.to.names(reference)
	}

	list.of.nonempty.spaces<-mixedsort(union(as.character(query@seqnames),as.character(reference@seqnames)))
	
	#to get here, chr is to have at least one represetative anywhere
	list.of.spaces<-intersect(list.of.nonempty.spaces,chromosomes.to.proceed)


	# now, list.of.spaces is what we will really proceed	
	if (length(intersect(as.character(query@seqinfo@seqnames),as.character(reference@seqinfo@seqnames)))==0)
		stop("There is no intersection between input set of common\n  chromosome names of query and reference and chromosomes.to.proceed.\n  It's all lost!")

	#if the awhole include list is empty, it is set to whole list
	if(length(chromosomes.to.include.in.awhole)==0) chromosomes.to.include.in.awhole<-list.of.spaces

	#we test whether all the chromosomes.to.include.in.awhole is in chromosomes.to.proceed	
	wrong_diff<-setdiff(chromosomes.to.include.in.awhole,list.of.spaces)
	if (length(wrong_diff)>0)
		stop("Some spaces  are in chromosomes.to.include.in.awhole but they are not in the process list.\n  It's all lost!\n",paste(wrong_diff,collapse="\n"), "\n")

	awhole.chromosomes<-intersect(chromosomes.to.include.in.awhole,list.of.spaces)

	wrong_diff<-setdiff(chromosomes.to.exclude.from.awhole,chromosomes.to.include.in.awhole)
	if (length(wrong_diff)>0)
		stop("Some spaces are in chromosomes.to.exclude.from.awhole but they are not in the awhole inclusion list.\n  It's all lost!:\n",paste(wrong_diff,collapse="\n"), "\n")
	
	#alternative is to be initial substring of 'attraction', 'repulsion' or 'two.sided'
	alt<-NA
	for (test in c("attraction","repulsion","two.sided")) {
		if (str_detect(test,str_c("^",alternative))) {alt<-test}
	}

	if(is.na(alt)) {stop("Alternative is not 'attraction' or 'repulsion' or 'two.sided'. \n  It's all lost!" )}

	#now, everyting is cheched
	awhole.chromosomes<-setdiff(chromosomes.to.include.in.awhole,chromosomes.to.exclude.from.awhole)

	if (length(awhole.chromosomes) < 2) 
	#cannot imagine what for awhole.chromosomes can have only one chromosome
	{
		awhole.chromosomes<-c()
#		warning(paste("Your awhole.chromosomes list contains only one chromosomes",awhole.chromosomes,"\n  so I switch the awhole information gather off."))
	}

	if (awhole.only && length(awhole.chromosomes) == 0)
		stop("The awhole.chromosomes list is empty and awhole.only switch is on.\n  It's all lost!")

	# we do not intereasted in lengths of those chromosomes that are
	# not in intersection of spacesA and spacesB;
	# if we have information in chromosomes.length, we do not care
	# about what was given in seqlengths of A and B
	result<-.GRangesGenometricsCorrelation(
			query=query,
			reference=reference,
			list.of.spaces=list.of.spaces,
			map.to.half=map.to.half,
			showProgressBar=showProgressBar,
			showTkProgressBar=showTkProgressBar,
			chromosomes.length=chromosomes.length,
			suppress.evaluated.length.warning=suppress.evaluated.length.warning,
			ecdf.area.permut.number=ecdf.area.permut.number,
			mean.distance.permut.number=mean.distance.permut.number,
			jaccard.measure.permut.number=jaccard.measure.permut.number,
			jaccard.permut.is.rearrangement=jaccard.permut.is.rearrangement,
			alternative=alt, #alt is the result of the check of alternative
			awhole.chromosomes=awhole.chromosomes,
			awhole.space.name=awhole.space.name,
			awhole.only=awhole.only,
			keep.distributions=keep.distributions,
			query.representing.point.function=query.representing.point.function,
			reference.representing.point.function=reference.representing.point.function
	)

	result@config$data=list()
	result@config$data$query=''
	result@config$data$reference=''
	result@config$chromosomes=vector('list',length(list.of.spaces)) 
	# empty list; vector('list') and list() is the same
	names(result@config$chromosomes)=list.of.spaces
	result@config$options=list()
	result@config$options$add.chr.as.prefix=add.chr.as.prefix
	result@config$options$awhole.only=awhole.only
	result@config$options$showProgressBar=showProgressBar
	result@config$options$showTkProgressBar=showTkProgressBar
	result@config$chromosomes.length=as.list(chromosomes.length)
	result@config$options$suppress.evaluated.length.warning=suppress.evaluated.length.warning
	result@config$options$keep.distributions=keep.distributions
	if (awhole.space.name!="awhole") result@config$options$awhole.space.name=awhole.space.name
	result@config$tests=list()
	result@config$tests$permut.number=permut.number
	result@config$tests$alternative=alt #alt is the result of the check of alternative
	if (ecdf.area.permut.number!=permut.number) result@config$tests$ecdf.area.permut.number=ecdf.area.permut.number
	if (mean.distance.permut.number!=permut.number) result@config$tests$mean.distance.permut.number=mean.distance.permut.number
	if (jaccard.measure.permut.number!=permut.number) result@config$tests$jaccard.measure.permut.number=jaccard.measure.permut.number
	#print(result@config)

	return(result)
}


.space_ranges<-function(granges,seqname){(plyranges::filter(granges, seqnames==seqname))@ranges}

.GRangesGenometricsCorrelation<-function(
	query,reference,
	list.of.spaces,
	map.to.half=TRUE,
	showProgressBar=TRUE,
	showTkProgressBar=FALSE,
	chromosomes.length=c(),
	suppress.evaluated.length.warning=FALSE,
	ecdf.area.permut.number=100,
	mean.distance.permut.number=100,
	jaccard.measure.permut.number=100,
	jaccard.permut.is.rearrangement=FALSE,
	alternative="two.sided",
	awhole.chromosomes=c(),
	awhole.space.name="awhole",
	awhole.only=FALSE,
	keep.distributions=FALSE,
	query.representing.point.function,
	reference.representing.point.function
)
{
	#the thing actually calculates everything
	#it is to called from GenomertiCorrelation
	#query,reference are two GRanges that are under analysis
	#list.of.spaces is list of spaces to work with
	#map.to.half is whether to calculate relative distances in [0,0.5] or in [0,1]
	#showProgressBar whether to show a progress indicator
	#showTkProgressBar whether to show a Tk progress indicator; work only if tcltk is loaded
	#chromosomes.length is an array of length of chromosomes, names are the names of chromosomes (we use it historically)
	#ecdf.area.permut.number is number of permutations for ecdf area method
	#mean.distance.permut.number is the same thing about the mean ref-to-query distance
	#awhole.chromosomes is list of chromosomes to be included in whole-genome calculation
	#awhole.space.name is the space name for this operation
	
	if (length(awhole.chromosomes) > 0)
		do.awhole<-TRUE
	else 
		do.awhole<-FALSE
	
	# HJ -- tcltk now part of R; we can Depends it
	if ((showTkProgressBar) && !require("tcltk",quietly=TRUE))
		showTkProgressBar=FALSE

	#if there is no loadable tcltk, switch showTkProgressBar off 

	if (showProgressBar || showTkProgressBar)
	{
		
		pb_capacity<-
			length(list.of.spaces)*8+ #in chr steps; non-permut
			length(list.of.spaces)*(ecdf.area.permut.number+mean.distance.permut.number+jaccard.measure.permut.number)+ #permutation
			3 #gather results
		if (do.awhole)
			pb_capacity<-pb_capacity+1+ecdf.area.permut.number+mean.distance.permut.number+jaccard.measure.permut.number
	}
	#indicator step : if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(pb)[1]+1)
	
	if (showProgressBar)
		txt_pb <- txtProgressBar(min = 0, max = pb_capacity,initial=0, style = 3)
		# create progress bar
	#indicator step : if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)


	if (showTkProgressBar)
	{
		tk_pb <- tkProgressBar(min = 0, max = pb_capacity,initial=0)
		done<-0;
		done_info<-'Starting...'
		#sprintf('%i of %i done',done,pb_capacity)
		setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		# create progress bar
	}
	#indicator step : 
	#if (showTkProgressBar)
	#{
	#	done <- getTkProgressBar(tk_pb)[1]+1;
	#	done_info <- sprintf('%i of %i done',done,pb_capacity)
	#	setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
	#}
	
	#we need to know all the chrom lengths or at least to mark it as NA
	for ( space in list.of.spaces )
	{
		my_space_length<-NA
		if (!is.na(seqlengths(query)[space])){
			my_space_length=seqlengths(query)[space]
		}
		if (!is.na(seqlengths(reference)[space])){
			if (!is.na(my_space_length) && my_space_length != seqlengths(reference)[space]) {
				stop(sprintf("Different length in query and reference for chromosone %s.\n",space),call.=FALSE)
			}
			my_space_length=seqlengths(reference)[space]
		}
		
		if (space %in% names(chromosomes.length) && !is.na(chromosomes.length[space]))
			my_space_length=chromosomes.length[space]
	
		if ( is.na(my_space_length) )
		{
			que_ranges<-.space_ranges(query,space)
			ref_ranges<-.space_ranges(reference, space)
			my_space_length<-chromosomes.length.eval(que_ranges,ref_ranges)
			if (! (suppress.evaluated.length.warning))
				warning(paste0("Length for chromosome ",space," is evaluated as ",as.character(chromosomes.length[space])," rather than pre-given."))
		}
		chromosomes.length[space]=my_space_length
		#here, we tested ends not to stick put, but now -- it is GRanges, it is already tested

	}	
	
	#initialise it all
	result<-list()
	
	
	#starting to prepare result
	for (space in list.of.spaces)
	{
		result[[space]]<-list()
	}

	if (do.awhole)
	{
		result[[awhole.space.name]]<-list()

		result[[awhole.space.name]][['query.population']]<-0
		result[[awhole.space.name]][['reference.population']]<-0
		result[[awhole.space.name]][['alternative']]<-alternative
		result[[awhole.space.name]][['query.coverage']]<-0
		result[[awhole.space.name]][['reference.coverage']]<-0
		result[[awhole.space.name]][['projection.test']]<-c()
		result[[awhole.space.name]][['projection.test']][['space.length']]<-0
		result[[awhole.space.name]][['projection.test']][['reference.coverage']]<-0
		result[[awhole.space.name]][['projection.test']][['query.hits']]<-0
		result[[awhole.space.name]][['query.reference.intersection']]<-0
		result[[awhole.space.name]][['query.reference.union']]<-0
		result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]<-0
		#makes no sense if no (mean.distance.permut==0) relative distance permutations done
		if (mean.distance.permut.number) 
		{
			#we initialise it here because we are goin to calculate it as a chomosomewise sum
			result[[awhole.space.name]][['scaled.absolute.min.distance.sum.null.list']]<-c()
		}
		# we initialise relative.distances.ecdf.deviation.area.null.list 
		# and jaccard.measure.null.list and jaccard.intersection.null.list later, 
		# in the do_awhole after main chromwise cycle 
		result[[awhole.space.name]][['relative.distances.data']]<-c()

		if (keep.distributions) # we need it only if we want distribution
		{
			result[[awhole.space.name]][['absolute.min.distance.data']]<-c()
	 		result[[awhole.space.name]][['absolute.inter.reference.distance.data']]<-c()
		}
	}
	else
		awhole.space.name<-c()
	
	
	#Kolmogorov-Smirnov test (local, pointset-to-pointset)
	#need: list of relative distances
	#combine to awhole: combine list
	#corr_sign is the same

	#Bernoulli test (global, pointset-to-ranges)	
	#need: ref_chromosomes.length; coverage of reference; 
	#number of query; number of query 
	#point fall into refrence
	#genomewide - just summarize these 4 chromosomewide

	#Relative distances between ECDF and diagonal (local, pointset-to-pointset)
	#p-value by permutations, permuting just a list of values in [0,0.5] of the same length as number of the query
	#genomewide - permute chromosome after chromosome, combine it into genomewide permuted ECDF
	#calculate the permuted area for each chromosome and for the genome awhole
	

	#Mean absolute query-to-reference distance (global, pointset-to-pointset)
	#p-value by permutation, the thing to be permuted is a set of positions of query points.
	#then, the mean value of query-ref absolute distance is calculated
	#we scale each chromosome's by length/#ref points
	#genomewide - permute chromosome after chromosome, 
	#then calculate the common mean of all the data - it is the genomewide value
	#then, the mean for each chromosome and for the genome yeilds a p-value

	#jaccard=query.reference.intersection/query.reference.union
	#calculate for each and for awhole
	#permute: each chrom and awhole is a sum for each perm (like all permutaqtions)
	#permutation can be permutation or rearrangement (user choice)


	if (map.to.half) rel.dist.top=0.5 else rel.dist.top=1.;

	for (space in list.of.spaces) #calculate everything nonpermutted for each chromosomes
	{
		qu<-sorted.representing.points(
			iranges=.space_ranges(query,space),
			representing.point.function=query.representing.point.function,
			chromosome.length=chromosomes.length[space],
			space=space
		)
		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}

		ref<-sorted.representing.points(
			iranges=.space_ranges(reference,space),
			representing.point.function=reference.representing.point.function,
			chromosome.length=chromosomes.length[space],
			space=space
		)

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}

		result[[space]][['query.population']]<-length(qu)

		result[[space]][['reference.population']]<-length(ref)

		result[[space]][['alternative']]<-alternative

		result[[space]][['query.coverage']]<-sum(width(reduce(.space_ranges(query,space))))

		result[[space]][['reference.coverage']]<-sum(width(reduce(.space_ranges(reference,space))))

		qu_or_ref_is_empty<-(length(qu)==0 || length(ref)==0)
		
		if (!qu_or_ref_is_empty) {
			result[[space]][['relative.distances.data']]<-
				query_to_ref_relative_distances(
						qu,ref,map.to.half,
						is_query_sorted=T,
						is_ref_sorted=T,
						chrom_length=chromosomes.length[space]
				)
		} else {
			result[[space]][['relative.distances.data']]<-c()
		}
		
		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
		if (!qu_or_ref_is_empty) {
			result[[space]][['relative.distances.ks.p.value']]<-
				ks.test(untie(result[[space]]$relative.distances.data),
					punif,min=0,max=rel.dist.top)$p.value
		} else {
			result[[space]][['relative.distances.ks.p.value']] <- 1.
		} #empty distribution

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
		
		if (!qu_or_ref_is_empty) {
			result[[space]][['relative.distances.ecdf.deviation.area']]<-
				integrate(
					function(x){
						return(abs((ecdf(result[[space]]$relative.distances.data))(x)-
						x/rel.dist.top))
					},
					lower=0,upper=rel.dist.top,
					subdivisions=length(result[[space]]$relative.distances.data)*100,
					rel.tol=integr_rel_tol)$value
		} else {
			result[[space]][['relative.distances.ecdf.deviation.area']]<-0
		}
		#makes no sense if no (ecdf.area.permut.number==0) relative distance permutations done
		if (ecdf.area.permut.number>0) 
		{
			result[[space]][['relative.distances.ecdf.deviation.area.null.list']]<-c()
		}

		relative.ecdf.area<-rel.dist.top-mean(result[[space]]$relative.distances.data)
			#rel.dist.top is integral(1) over (0..rel.dist.top)
			
			#integrate(
			#		ecdf(result[[space]]$relative.distances.data),
			#		lower=0,upper=rel.dist.top,
			#		subdivisions=length(result[[space]]$relative.distances.data)*100,
			#		rel.tol=integr_rel_tol)$value
		indep_area=rel.dist.top/2
		result[[space]][['relative.distances.ecdf.area.correlation']]=(relative.ecdf.area-indep_area)/indep_area
		#it is Exp(dist)

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}


		result[[space]][['query.reference.intersection']]<-sum(width(reduce(intersect(.space_ranges(query,space),.space_ranges(reference,space)))))
		result[[space]][['query.reference.union']]<-sum(width(reduce(union(.space_ranges(query,space),.space_ranges(reference,space)))))

		if (result[[space]][['query.reference.union']] > 0)
			result[[space]][['jaccard.measure']]<-result[[space]][['query.reference.intersection']]/result[[space]][['query.reference.union']]	
		else
			result[[space]][['jaccard.measure']]<-0

		if (jaccard.measure.permut.number>0) 
		{
			result[[space]][['jaccard.measure.null.list']]<-c()
			result[[space]][['jaccard.intersection.null.list']]<-c()
		}
		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
		
		result[[space]][['projection.test']]<-
				  query.to.ref.projection.statistics(qu,.space_ranges(reference,space),TRUE,chromosomes.length[space])	
				  
		if(!qu_or_ref_is_empty) {
			result[[space]][['projection.test.obs.to.exp']]<- 
				(	
					result[[space]][['projection.test']][['query.hits']]
					/
					result[[space]][['query.population']]
				)*
				(
					result[[space]][['projection.test']][['space.length']] 
					/ 
					result[[space]][['projection.test']][['reference.coverage']]
				)
			
			if(alternative == "two.sided") {
				if (result[[space]][['projection.test.obs.to.exp']] < 1.) {#repulsion
						direction<-"repulsion"
				} else {
						direction<-"attraction"
				}
			}	
			else {direction<-alternative}

			proj.p.value<-
				pbinom(
					result[[space]][['projection.test']][['query.hits']],
					result[[space]][['query.population']],
					result[[space]][['projection.test']][['reference.coverage']]/
								result[[space]][['projection.test']][['space.length']],
					lower.tail = direction == "repulsion")
			
			if(alternative == "two.sided") {
				 proj.p.value<- min(proj.p.value*2,1.)
				 result[[space]][['projection.test.direction']]<-direction
			}
			result[[space]][['projection.test.p.value']]<-proj.p.value
		} else {
			result[[space]][['projection.test.p.value']] <- 1.
			result[[space]][['projection.test.direction']] <- "undefined"
			result[[space]][['projection.test.obs.to.exp']] <- 1.
		}
		
		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}

		result[[space]][['absolute.min.distance.data']]<-
			query_to_ref_min_absolute_distances(
				qu,ref,map.to.half,is_query_sorted=T,is_ref_sorted=T,chrom_length=chromosomes.length[space]
		)
		result[[space]][['scaled.absolute.min.distance.sum']]<-
			sum(result[[space]][['absolute.min.distance.data']])*length(ref)/(chromosomes.length[[space]])

		if(mean.distance.permut.number>0 || keep.distributions) 
		# makes sense only if permutation test on absolute distances test is ordered or if we want to keep
		# the distribution of absolute min distances 
		{
			result[[space]][['absolute.inter.reference.distance.data']]<-ref[2:length(ref)]-ref[1:length(ref)-1]
			#collect all ref-ref distances
			result[[space]][['absolute.inter.reference.distance.data']]<-c(result[[space]][['absolute.inter.reference.distance.data']],chromosomes.length[[space]]+ref[1]-ref[length(ref)])	
			#add the circilar one	

			#makes no sense if no (mean.distance.permut==0) absolute distance permutations done
			
			if (mean.distance.permut.number>0) 
			{
				result[[space]][['reference.middles']]<-ref #we will need it for permutations 
					
				result[[space]][['scaled.absolute.min.distance.sum.null.list']]<-c()
				#result[[space]][['unscaled.absolute.min.distance.sum.null.list']]<-c()
			}
		}

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
		#cat('\n#',space,'#\n')

		#put it to the whole-genome distribution
		if (do.awhole && space %in% awhole.chromosomes)
		{
			result[[awhole.space.name]][['relative.distances.data']]<-
						c(result[[awhole.space.name]][['relative.distances.data']],result[[space]][['relative.distances.data']])

			result[[awhole.space.name]][['absolute.min.distance.data']]<-
						c(result[[awhole.space.name]][['absolute.min.distance.data']],result[[space]][['absolute.min.distance.data']])

			result[[awhole.space.name]][['absolute.inter.reference.distance.data']]<-
						c(result[[awhole.space.name]][['absolute.inter.reference.distance.data']],result[[space]][['absolute.inter.reference.distance.data']])

			for (stats in names(result[[space]][['projection.test']]))
					result[[awhole.space.name]][['projection.test']][[stats]]<-
							result[[awhole.space.name]][['projection.test']][[stats]]+result[[space]][['projection.test']][[stats]]


			result[[awhole.space.name]][['query.reference.intersection']]<-result[[awhole.space.name]][['query.reference.intersection']]+result[[space]][['query.reference.intersection']]
			result[[awhole.space.name]][['query.reference.union']]<-result[[awhole.space.name]][['query.reference.union']]+result[[space]][['query.reference.union']]

			result[[awhole.space.name]][['query.population']]<-result[[awhole.space.name]][['query.population']]+result[[space]][['query.population']]
			result[[awhole.space.name]][['reference.population']]<-result[[awhole.space.name]][['reference.population']]+result[[space]][['reference.population']]

			result[[awhole.space.name]][['query.coverage']]<-result[[awhole.space.name]][['query.coverage']]+result[[space]][['query.coverage']]
			result[[awhole.space.name]][['reference.coverage']]<-result[[awhole.space.name]][['reference.coverage']]+result[[space]][['reference.coverage']]
			#makes sense even if no (mean.distance.permut==0) absolute distance permutations done
			#if(mean.distance.permut.number>0) 
			#{
				result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]<-
					result[[awhole.space.name]][['scaled.absolute.min.distance.sum']]+
					result[[space]][['scaled.absolute.min.distance.sum']]
			#}
		}
	}
	
	#whole genome results started
	#relative distances data for whole genome
	if (do.awhole)
	{
		space<-awhole.space.name	
		the.names<- names(result[[space]])

		#qu_or_ref_is_empty<-(result[[space]][['query.population']]==0 ||result[[space]][['reference.population']]==0 )

		if ('relative.distances.data' %in% the.names)
		{
			result[[space]][['relative.distances.ks.p.value']]<-
				ks.test(untie(result[[space]]$relative.distances.data),punif,min=0,max=rel.dist.top)$p.value

			result[[space]][['relative.distances.ecdf.deviation.area']]<-
				integrate(
						function(x){return(abs((ecdf(result[[space]]$relative.distances.data))(x)-x/rel.dist.top))},
						lower=0,upper=rel.dist.top,
						subdivisions=length(result[[space]]$relative.distances.data)*100,
						rel.tol=integr_rel_tol)$value
			
			#makes no sense if no (ecdf.area.permut.number==0) permutations done
			if (ecdf.area.permut.number>0) 
			{
				result[[space]][['relative.distances.ecdf.deviation.area.null.list']]<-c()
			}

			relative.ecdf.area<-rel.dist.top-mean(result[[space]]$relative.distances.data)
				#rel.dist.top is integral(1) over (0..rel.dist.top)
				
				#integrate(
				#		ecdf(result[[space]]$relative.distances.data),
				#		lower=0,upper=rel.dist.top,
				#		subdivisions=length(result[[space]]$relative.distances.data)*100,
				#		rel.tol=integr_rel_tol)$value
			indep_area=rel.dist.top/2

			result[[space]][['relative.distances.ecdf.area.correlation']]=(relative.ecdf.area-indep_area)/indep_area

		}

		
		if ('projection.test' %in% the.names) #projection test for awhole genome
		{
			if(!qu_or_ref_is_empty) {
				result[[space]][['projection.test.obs.to.exp']]<- 
					(	
						result[[space]][['projection.test']][['query.hits']]
						/
						result[[space]][['query.population']]
					)*
					(
						result[[space]][['projection.test']][['space.length']] 
						/ 
						result[[space]][['projection.test']][['reference.coverage']]
					)
				
				if(alternative == "two.sided") {
					if (result[[space]][['projection.test.obs.to.exp']] < 1.) {#repulsion
							direction<-"repulsion"
					} else {
							direction<-"attraction"
					}
				}	
				else {direction<-alternative}

				proj.p.value<-
					pbinom(
						result[[space]][['projection.test']][['query.hits']],
						result[[space]][['query.population']],
						result[[space]][['projection.test']][['reference.coverage']]/
									result[[space]][['projection.test']][['space.length']],
						lower.tail = direction == "repulsion")
				
				if(alternative == "two.sided") {
					 proj.p.value<- min(proj.p.value*2,1.)
					 result[[space]][['projection.test.direction']]<-direction
				}
				result[[space]][['projection.test.p.value']]<-proj.p.value
			} else {
				result[[space]][['projection.test.p.value']] <- 1.
				result[[space]][['projection.test.direction']] <- "undefined"
				result[[space]][['projection.test.obs.to.exp']] <- 1.
			}
		}

		if ('query.reference.union' %in% the.names && 'query.reference.intersection' %in% the.names)
		{
			if (result[[space]][['query.reference.union']] > 0)
				result[[space]][['jaccard.measure']]<-result[[space]][['query.reference.intersection']]/result[[space]][['query.reference.union']]	
			else
				result[[space]][['jaccard.measure']]<-0

			if (jaccard.measure.permut.number>0)
			{
				result[[space]][['jaccard.measure.null.list']]<-c()
				result[[space]][['jaccard.intersection.null.list']]<-c()
			}

		}

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',awhole.space.name,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
	}

	#ecdf permutation gathering
	if (ecdf.area.permut.number>0)
		for (permutaion.number in 1:ecdf.area.permut.number)
		{
			awhole_sample<-c()
			for (space in list.of.spaces)
			{
				#each choromosome
				sample_size<-length(.space_ranges(query,space))
				sample<-runif(sample_size,0,rel.dist.top) #sample
				if (length(sample)>0) {
				dev_area<-
					integrate(
							function(x){return(abs((ecdf(sample)(x))-x/rel.dist.top))},
							lower=0,upper=rel.dist.top,
							subdivisions=sample_size*100,
							rel.tol=integr_rel_tol)$value
				} else {dev_area<-0}
				#so if query is empty it returns 0
				#if the reference is - we do not care here

				result[[space]][['relative.distances.ecdf.deviation.area.null.list']]<-
					c(result[[space]][['relative.distances.ecdf.deviation.area.null.list']],dev_area) #save the result

				if (space %in% awhole.chromosomes)
					awhole_sample<-c(awhole_sample,sample) # put the sample to awole_genome_sample
				if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
				if (showTkProgressBar)
				{
					done <- getTkProgressBar(tk_pb)[1]+1;
					done_info <- sprintf('ecdf permut. #%i ; %i of %i done',permutaion.number,done,pb_capacity)
					setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
				}
			}
			#now, this permutaion for the whole genome
			if (!do.awhole) next
			space<-awhole.space.name
			dev_area<-
				integrate(
						function(x){return(abs((ecdf(awhole_sample)(x))-x/rel.dist.top))},
						lower=0,upper=rel.dist.top,
						subdivisions=length(awhole_sample)*100,
						rel.tol=integr_rel_tol)$value
			result[[space]][['relative.distances.ecdf.deviation.area.null.list']]<-
				c(result[[space]][['relative.distances.ecdf.deviation.area.null.list']],dev_area) #save the result
			if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
			if (showTkProgressBar)
			{
				done <- getTkProgressBar(tk_pb)[1]+1;
				done_info <- sprintf('ecdf permut. #%i ; %i of %i done',permutaion.number,done,pb_capacity)
				setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
			}
		}

	#mean min distance permutation gathering
	if (mean.distance.permut.number>0 )
	{
		for (permutaion.number in 1:mean.distance.permut.number)
		{
			awhole_scaled_sum<-0
			for (space in list.of.spaces)
			{
				#each choromosome
				sample_size<-length(.space_ranges(query,space))
				chr_length<-chromosomes.length[[space]]
				ref<-result[[space]][['reference.middles']]
				sample<-((floor(runif(sample_size,0,1)*chr_length*2)+1)/2) #sample
				mean_min_sum<-query_to_ref_min_absolute_distance_sum(
					sample,ref,
					map.to.half,
					is_query_sorted=FALSE,
					is_ref_sorted=TRUE,
				chrom_length=chr_length)
				scaled_mean_min_sum<-mean_min_sum*(length(ref)/chr_length)
				result[[space]][['scaled.absolute.min.distance.sum.null.list']]<-
					c(result[[space]][['scaled.absolute.min.distance.sum.null.list']],scaled_mean_min_sum) #save the result

				#result[[space]][['unscaled.absolute.min.distance.sum.null.list']]<-
				#	c(result[[space]][['unscaled.absolute.min.distance.sum.null.list']],mean_min_sum) #save the result
				if (space %in% awhole.chromosomes)
					awhole_scaled_sum<-awhole_scaled_sum+scaled_mean_min_sum # put the sample to awole_genome_sample
				if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
				if (showTkProgressBar)
				{
					done <- getTkProgressBar(tk_pb)[1]+1;
					done_info <- sprintf('absdist permut. #%i; %i of %i done',permutaion.number,done,pb_capacity)
					setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
				}
			}
			#now, this permutaion for the whole genome
			if (!do.awhole) next
			space<-awhole.space.name
			result[[space]][['scaled.absolute.min.distance.sum.null.list']]<-
				c(result[[space]][['scaled.absolute.min.distance.sum.null.list']],awhole_scaled_sum) #save the result
			if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
			if (showTkProgressBar)
			{
				done <- getTkProgressBar(tk_pb)[1]+1;
				done_info <- sprintf('absdist permut. #%i; %i of %i done',permutaion.number,done,pb_capacity)
				setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
			}
		}
		result[[space]][['reference.middles']]<-NULL
		#we do need it any more
	}
	#jaccard.measure permutation gathering
	if (jaccard.measure.permut.number>0 )
		for (permutaion.number in 1:jaccard.measure.permut.number)
		{
			awhole.union<-0
			awhole.intersection<-0
			for (space in list.of.spaces)
			{
				#each choromosome
				chr_length<-chromosomes.length[[space]]
				if (!jaccard.permut.is.rearrangement)
					permuted.query<-.permuteIRanges(.space_ranges(query,space),chr_length)
				else
					permuted.query<-.rearrangeIRanges(.space_ranges(query,space),chr_length)

				space.intersection<-sum(width(reduce(intersect(permuted.query,.space_ranges(reference,space)))))
				space.union<-sum(width(reduce(union(permuted.query,.space_ranges(reference,space)))))
				
				if (space.union> 0)
					space.jaccard.measure<-space.intersection/space.union	
				else
					space.jaccard.measure<-0

				result[[space]][['jaccard.measure.null.list']]<-
					c(result[[space]][['jaccard.measure.null.list']],space.jaccard.measure)

				result[[space]][['jaccard.intersection.null.list']]<-
					c(result[[space]][['jaccard.intersection.null.list']],space.intersection) 
				
				if (space %in% awhole.chromosomes)
				{
					awhole.union<-awhole.union+space.union
					awhole.intersection<-awhole.intersection+space.intersection
				}
				if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
				if (showTkProgressBar)
				{
					done <- getTkProgressBar(tk_pb)[1]+1;
					done_info <- sprintf('Jaccard permut. #%i; %i of %i done',permutaion.number,done,pb_capacity)
					setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
				}
			}
			#now, this permutaion for the whole genome
			if (!do.awhole) next

			if (awhole.union> 0)
				awhole.jaccard.measure<-awhole.intersection/awhole.union	
			else
				awhole.jaccard.measure<-0
			space<-awhole.space.name
			result[[space]][['jaccard.measure.null.list']]<-
				c(result[[space]][['jaccard.measure.null.list']],awhole.jaccard.measure)
			result[[space]][['jaccard.intersection.null.list']]<-
				c(result[[space]][['jaccard.intersection.null.list']],awhole.intersection) 
			if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
			if (showTkProgressBar)
			{
				done <- getTkProgressBar(tk_pb)[1]+1;
				done_info <- sprintf('Jaccard permut. #%i; %i of %i done',permutaion.number,done,pb_capacity)
				setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
			}
		}
	if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
	if (showTkProgressBar)
	{
		done <- getTkProgressBar(tk_pb)[1]+1;
		done_info <- sprintf('gathering results; %i of %i done',done,pb_capacity)
		setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
	}


	#area p-values
	if (ecdf.area.permut.number>0)
	{
		for (space in c(list.of.spaces,awhole.space.name))
		{
			#if ( result[[space]][['query.population']]==0 || result[[space]][['reference.population']]==0) {
			#	p.value<-1.
			#} else {
			#	p.value<-
			#		1-(ecdf(result[[space]][['relative.distances.ecdf.deviation.area.null.list']]))(
			#			result[[space]][['relative.distances.ecdf.deviation.area']]
			#		) #we treat only right side
			#	if (p.value==0) 
			#		p.value = paste("<",toString(1/ecdf.area.permut.number),sep='')
			#}
			#result[[space]][['relative.distances.ecdf.deviation.area.p.value']]<-p.value
			# old end	-- er are wrong with right side???? ecdf. think 
			#we will remove it if the following works


			if ( result[[space]][['query.population']]==0 || result[[space]][['reference.population']]==0) {
				p.value<-1.
				direction<-"undefined"	
			} else {
				p.value<-
					1 - (ecdf(result[[space]][['relative.distances.ecdf.deviation.area.null.list']]))(
						result[[space]][['relative.distances.ecdf.deviation.area']]) #we treat only right side
					
				if('attraction'==alternative) { #it is for lower tail of distance 
					direction<-'attraction'
				} else if ('repulsion'==alternative) { #upper tail
					direction <- 'repulsion'
					p.value <- 1-p.value
				} else { #two.sided
				if (p.value<0.5) {
						p.value <- p.value*2
						direction <- 'attraction'
					} else {
						p.value <- (1-p.value)*2
						direction <- 'repulsion'
					}
				}
				if (p.value==0) {
					if (alternative=='two.sided') {
						p.value = paste("<",toString(2/mean.distance.permut.number),sep='')
					} else {
						p.value = paste("<",toString(1/mean.distance.permut.number),sep='')
					}
				}
			}
			result[[space]][['relative.distances.ecdf.deviation.area.p.value']]<-p.value
			if('two.sided'==alternative) {result[[space]][['relative.distances.ecdf.deviation.area.test.direction']] <- direction}
		}
	}
	if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
	if (showTkProgressBar)
	{
		done <- getTkProgressBar(tk_pb)[1]+1;
		done_info <- sprintf('gathering results; %i of %i done',done,pb_capacity)
		setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
	}
	
	
	#min distance sum p-values
	if (mean.distance.permut.number>0)
	{
		for (space in c(list.of.spaces,awhole.space.name))
		{
			if ( result[[space]][['query.population']]==0 || result[[space]][['reference.population']]==0) {
				p.value<-1.
				direction<-"undefined"	
			} else {
				p.value<-
					(ecdf(result[[space]][['scaled.absolute.min.distance.sum.null.list']]))(result[[space]][['scaled.absolute.min.distance.sum']]) #we treat only right side
					
				if('attraction'==alternative) { #it is for lower tail of distance 
					direction<-'attraction'
				} else if ('repulsion'==alternative) { #upper tail
					direction <- 'repulsion'
					p.value <- 1-p.value
				} else { #two.sided
					if (p.value<0.5) {
						p.value <- p.value*2
						direction <- 'attraction'
					} else {
						p.value <- (1-p.value)*2
						direction <- 'repulsion'
					}
				}
				if (p.value==0) {
					if (alternative=='two.sided') {
						p.value = paste("<",toString(2/mean.distance.permut.number),sep='')
					} else {
						p.value = paste("<",toString(1/mean.distance.permut.number),sep='')
					}
				}
			}
			result[[space]][['scaled.absolute.min.distance.sum.p.value']]<-p.value
			if('two.sided'==alternative) {result[[space]][['scaled.absolute.min.distance.sum.test.direction']] <- direction} 
		}
	}
	
	#jaccard distance sum p-values
	if (jaccard.measure.permut.number>0)
	{
		for (space in c(list.of.spaces,awhole.space.name))
		{
			if ( result[[space]][['query.population']]==0 || result[[space]][['reference.population']]==0) {
				p.value<-1.
				direction<-"undefined"	
			} else {
				p.value<-
					(ecdf(result[[space]][['jaccard.measure.null.list']]))(result[[space]][['jaccard.measure']]) #we treat only right side
				
				if('attraction'==alternative) { #it is for upper tail of distance 
					direction<-'attraction'
					p.value <- 1-p.value
				} else if ('repulsion'==alternative) { #lower tail
					direction <- 'repulsion'
				} else { #two.sided
					if (p.value<0.5) {
						p.value <- p.value*2
						direction <- 'repulsion'
					} else {
						p.value <- (1-p.value)*2
						direction <- 'attraction'
					}
				}
				if (p.value==0) {
					if (alternative=='two.sided') {
						p.value = paste("<",toString(2/mean.distance.permut.number),sep='')
					} else {
						p.value = paste("<",toString(1/mean.distance.permut.number),sep='')
					}
				}
			}
			result[[space]][['jaccard.measure.p.value']]<-p.value
			if('two.sided'==alternative) {result[[space]][['jaccard.measure.test.direction']] <- direction }
		}
	}
	if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)
	if (showTkProgressBar)
	{
		done <- getTkProgressBar(tk_pb)[1]+1;
		done_info <- sprintf('gathering results; %i of %i done',done,pb_capacity)
		setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
	}
	#clear all
	if(!keep.distributions)
	{
		for (space in c(list.of.spaces,awhole.space.name))
		{
			result[[space]][['projection.test']]<-NULL
			result[[space]][['relative.distances.data']]<-NULL
			#result[[space]][['relative.distances.ecdf.deviation.area']]<-NULL
			result[[space]][['relative.distances.ecdf.deviation.area.null.list']]<-NULL
			result[[space]][['scaled.absolute.min.distance.sum.null.list']]<-NULL
			result[[space]][['absolute.min.distance.data']]<-NULL
			result[[space]][['absolute.inter.reference.distance.data']]<-NULL
			result[[space]][['jaccard.measure.null.list']]<-NULL
			result[[space]][['jaccard.intersection.null.list']]<-NULL
		}
	}
	
	#now, we want the order of records in result[[awhole.space.name]] to be the same as for all the chromosomes

	if (do.awhole)
	{
		result.awhole<-result[[awhole.space.name]];

		result[[awhole.space.name]]<-list();
		
		list.of.names.in.return<-names(result[[1]]);
		# 1 here is just bacause element 1 exist whatever;
		# what we want is just to have the same order of fields in
		# result[[awhole.space.name]] as in all chromosome results

		if (awhole.only)
		#if we want to have only the awhole result,
		#we are to clear result now and to prepare 
		#result[[awhole.space.name]] as a list()
		{
			result<-list();
			result[[awhole.space.name]]<-list();
		}

		for (name in list.of.names.in.return)
		{
			result[[awhole.space.name]][[name]]<-result.awhole[[name]]
		}
	}
	if (showProgressBar) 
	{ 
		close(txt_pb)
	}
	if (showTkProgressBar) 
	{ 
		close(tk_pb)
	}
	
	return(new('GenometriCorrResult',.Data=result))
}


sorted.representing.points<-function(iranges,representing.point.function,chromosome.length,space)
{
	#iranges is IRanges
	mids<-representing.point.function(start(iranges),end(iranges),chromosome.length,space)
	return(mids[order(mids)])
}

#input: Irange
#output: vector of represetating points


chromosomes.length.eval<-function(ir_query_space,ir_reference_space)
{
	#ir_query and ir_reference are IRange
	min_coord<-max(min(start(ir_query_space),start(ir_reference_space),end(ir_query_space),end(ir_reference_space)),0)
	#we prefer this way rather than obvious max
	#to equal the error if everything is located in say right and left telomeres
	#in other words, if everything starts at 10000 we will add 10000 to the max coord of any interval
	max_coord<-max(start(ir_query_space),start(ir_reference_space),end(ir_query_space),end(ir_reference_space))
	return(max_coord+min_coord)
}


query_to_ref_relative_distances<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA)
#calculate distances for a pair of positions (points)  vectors: query, ref
{
	if (length(query)==0 || length(ref)==0) {
		return (c())
	}
	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]
	rel_data<-c()

	#if chrom_length==0, we will loose all probes that fall out ref range, 
	#otherwise we provide the chrom length and so circle the reference
	refindex<-1
	firstref<-ref[1]
	lastref<-ref[length(ref)]
	for (quindex in 1:length(query))
	#we like to calcilate wheb current ref < current query <= next ref
	{
		qu<-query[quindex]
		if (qu <= firstref ) 
		{
			#we are still left from first ref 
			if (is.na(chrom_length)) next
			left_wing <- qu + chrom_length-lastref
			right_wing <- firstref-qu		
		}
		else if (lastref<qu)
		{
			#we are right from the last ref
			if (is.na(chrom_length)) next
			left_wing <- qu-lastref
			right_wing <- chrom_length-qu+firstref
		}
		else
		{
			#we are inside
			while (ref[refindex+1] < qu) refindex<-refindex+1
			#ref[refindex] < qu <= ref[refindex+1]
			left_wing <- qu-ref[refindex]
			right_wing <- ref[refindex+1]-qu
		}
		#here, we know both wings
		#resampling all borders
		if (left_wing==0 || right_wing==0) 
		{
			if (runif(1)>0.5)
			{
				right_wing<-right_wing+left_wing
				left_wing<-0
			}
			else
			{
				left_wing<-right_wing+left_wing
				right_wing<-0
			}
		}
		#cat("*** ",left_wing,'   ',right_wing,'\n')
		rel_distance<-left_wing/(left_wing+right_wing)

		#
		if (map.to.half)	rel_distance<-min(rel_distance,1-rel_distance)

		rel_data<-c(rel_data,rel_distance) #add measure
	}
	#so we add noise to the query
	return(rel_data)
} #query_to_ref_relative_distances<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=0)

query_to_ref_min_absolute_distances<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA)
#calculate distances for a pair of positions (points)  vectors: query, ref
{
	if (length(query)==0 || length(ref)==0) {
		return (c())
	}

	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]

	#if chrom_length==0, we will loose all probes that fall out ref range, 
	#otherwise we provide the chrom length and so circle the reference
	refindex<-1
	firstref<-ref[1]
	lastref<-ref[length(ref)]
	result<-c()
	
	for (quindex in 1:length(query))
	#we like to calcilate wheb current ref < current query <= next ref
	{
		qu<-query[quindex]
		if (qu <= firstref ) 
		{
			#we are still left from first ref 
			if (is.na(chrom_length)) next
			left_wing <- qu + chrom_length-lastref
			right_wing <- firstref-qu		
		}
		else if (lastref<qu)
		{
			#we are right from the last ref
			if (is.na(chrom_length)) next
			left_wing <- qu-lastref
			right_wing <- chrom_length-qu+firstref
		}
		else
		{
			#we are inside
			while (ref[refindex+1] < qu) refindex<-refindex+1
			#ref[refindex] < qu <= ref[refindex+1]
			left_wing <- qu-ref[refindex]
			right_wing <- ref[refindex+1]-qu
		}
		#here, we know both wings
		#resampling all borders
		if (left_wing==0 || right_wing==0) 
		{
			if (runif(1)>0.5)
			{
				right_wing<-right_wing+left_wing
				left_wing<-0
			}
			else
			{
				left_wing<-right_wing+left_wing
				right_wing<-0
			}
		}
		if (map.to.half)	abs_distance<-min(left_wing,right_wing)
		else abs_distance<-left_wing

		#
		result<-c(result,abs_distance) #add this
	}
	#so we add noise to the query
	return(result)
} #query_to_ref_average_mean_distance<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA) 

query_to_ref_min_absolute_distance_sum<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA)
#calculate distances for a pair of positions (points)  vectors: query, ref
{

	if (length(query)==0 || length(ref)==0) {
		return (0)
	}

	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]

	#if chrom_length==0, we will loose all probes that fall out ref range, 
	#otherwise we provide the chrom length and so circle the reference
	refindex<-1
	firstref<-ref[1]
	lastref<-ref[length(ref)]
	sum_of_min_distances<-0
	
	for (quindex in 1:length(query))
	#we like to calcilate wheb current ref < current query <= next ref
	{
		qu<-query[quindex]
		if (qu <= firstref ) 
		{
			#we are still left from first ref 
			if (is.na(chrom_length)) next
			left_wing <- qu + chrom_length-lastref
			right_wing <- firstref-qu		
		}
		else if (lastref<qu)
		{
			#we are right from the last ref
			if (is.na(chrom_length)) next
			left_wing <- qu-lastref
			right_wing <- chrom_length-qu+firstref
		}
		else
		{
			#we are inside
			while (ref[refindex+1] < qu) refindex<-refindex+1
			#ref[refindex] < qu <= ref[refindex+1]
			left_wing <- qu-ref[refindex]
			right_wing <- ref[refindex+1]-qu
		}
		#here, we know both wings
		#resampling all borders
		if (left_wing==0 || right_wing==0) 
		{
			if (runif(1)>0.5)
			{
				right_wing<-right_wing+left_wing
				left_wing<-0
			}
			else
			{
				left_wing<-right_wing+left_wing
				right_wing<-0
			}
		}
		if (map.to.half)	abs_distance<-min(left_wing,right_wing)
		else abs_distance<-left_wing

		#
		sum_of_min_distances<-sum_of_min_distances + abs_distance #add this
	}
	#so we add noise to the query
	return(sum_of_min_distances)
} #query_to_ref_average_mean_distance<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA) 

query.to.ref.projection.statistics<-function(query,ref,is_query_sorted=F,chrom_length=NA)
#query is vector of points; ref is an IRanges object
{
	if (is.na(chrom_length))
		stop("query.to.ref.projection.statistics cannnot work if chrom_length is NA")

	#if (!is_query_sorted) 
	#	query<-query[order(query)]

	#we reduce it, whatever!
	ref<-reduce(ref) #remove duplication in coverage

	#ref<-ref[order(ref)]
	
	projection_data<-c()
	projection_data[['space.length']]<-chrom_length
	projection_data[['reference.coverage']]<-sum(width(ref))
	
	#now, we revert our q to IRanges. If qu[i] in integer then the range has lentgh 1 and its start and end are q[i]
	#if q[i] is integer+0.5 we assign it to [integer,integer+1] interval

	ir_query<-IRanges(start=as.integer(query),width=ifelse(query %% 1 == 0,1,2))

	projection_data[['query.hits']]<-sum(overlapsAny(ir_query,ref))

	return(projection_data)
} #query.to.ref.projection.statistics<-function(query,ref,is_query_sorted=F,chrom_length=NA)


untie<-function(input_list)
{
	if (length(input_list)<=1) return (input_list)
	shaking_max <- epsilon
	input_list <- input_list[order(input_list)]
	for (i in (1:(length(input_list)-1)))
	{
		if (input_list[i]==input_list[i+1])
			input_list[i] <- input_list[i]+shaking_max*(-0.5+runif(1))
	}
	return(input_list)
}
#untie remove ties (i.e. equal values) from a data vector

