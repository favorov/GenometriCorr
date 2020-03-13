# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2020 Alexander Favorov, Loris Mularoni, 
#               Yulia Medvedeva, Harris A. Jaffee, 
#               Ekaterina V. Zhuravleva, Veronica Busa,
#               Leslie M. Cope, Andrey A. Mironov, 
#               Vsevolod J. Makeev, Sarah J. Wheelan.
#
# GenomertiCorrelation is the main function of the package


#' MapRangesToGenomicIntervals
#'
#' MapRangesToGenomicIntervals collapses all the chromosomes to their subsets covered by where.to.map and then liftover the what.to.map to this genome subset. The function is intended to prepare an estimation correlation of two features that are spatially restricted in the genome. We map them both to their location areas (the \code{where.to.map} is common for query and for the rerefence), and then we calculate the \code{GenometriCorrelation} between the mapped fatures.)
#'
#' @param where.to.map	The set of genomic intervals we map to.
#' @param what.to.map  The set of ranges that we map.
#' @param chrom.suffix The suffix to be appended to all the sestination chromosome names in the mapping; default is "_mapped".
#' @param chromosomes.to.proceed The default set of chromosomes to map is the intersection of the chromosomes in where.to.map and what.to.map. If we want to restrict the set, we can do it with this parameter.
#' @param chromosomes.length is an alternative to seqingo() ot the GRanges of where.to.map way to pass the lengths of chromosomes to the mapping routine
#' @param unmapped.chromosome.warning For each chromosome that is represented in \code{what.to.map} and that is included in \code{chromosomes.to.proceed} if it is given and that is not represented in \code{where.to.map}, a warning is generated if \code{unmapped.chromosome.warning} is \code{TRUE}. The default is \code{TRUE}.
#' @param nonnormalised.mapping.warning	If the input mapping space is not normalised (e.g. contains overlapping intervals), it is normalised before the mapping. A warning is generated if \code{nonnormalised.mapping.warning} is \code{TRUE}. The default is \code{TRUE}.
#' @return GRanges object that is a liftover of \code{what.to.map} to subgenome covered by \code{where.to.map}  
#' @author Alexander Favorov \email{favorov@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Veronica Busa, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#' @references \href{http://genometricorr.sourceforge.net/}{GenometriCorr home}
#' @seealso The \code{\link{GenometriCorr}} documentation and vignette.
#' @examples
#' library('GenometriCorr')
#' intervals<-GRanges(ranges=IRanges(c(1,10001,1,10001),width=c(1000)),seqnames=c('chr1','chr1','chr2','chr2'),seqlengths=c('chr1'=300000,'chr2'=300000))
#' ranges=GRanges(ranges=IRanges(c(10,110,10101,100000,500,550,1055),width=c(10)),seqnames=c(rep('chr1',4),rep('chr2',3)),seqlengths=c('chr1'=300000,'chr2'=300000))
#' mapped<-MapRangesToGenomicIntervals(where.to.map=intervals,what.to.map=ranges)
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

#' @keywords manip 
#' @export
#' @import rtracklayer plyranges
MapRangesToGenomicIntervals<-function(
	where.to.map, what.to.map,
	chrom.suffix="_mapped",
	chromosomes.to.proceed=NA,
	chromosomes.length=c(),
	unmapped.chromosome.warning=TRUE,
	nonnormalised.mapping.warning=TRUE
)
#subgenome and are suppose to be GRanges
#in the second case, we test the lenthg equivalence
{
	is.where.gr<-inherits(where.to.map,"GRanges")
	if (!is.where.gr)
		stop("where.to.map is not GRanges. It's all lost!\n")	
	chromosome.names.where<-NA
	chromosome.names.where<-as.character(unique(seqnames(where.to.map)))

	is.what.gr<-inherits(what.to.map,"GRanges")
	if (!is.what.gr)
		stop("what.to.map is not GRanges. It's all lost!\n")	
	
	#we are here, they are granges
	if (!is.na(chromosomes.to.proceed)) {
		where.to.map <- where.to.map %>% filter(seqnames %in% chromosomes.to.proceed)
		what.to.map <- what.to.map %>% filter(seqnames %in% chromosomes.to.proceed)
	}
	
	if(! identical(where.to.map,reduce(where.to.map))){
    if (nonnormalised.mapping.warning) {warning("GRanges object to map to includes overlapping intervals. Normalised.")}
		where.to.map<-reduce(where.to.map)
  }

	if(unmapped.chromosome.warning) {
		unmapped_chroms<-setdiff(what.to.map@seqnames,where.to.map@seqnames)
		if(length(unmapped_chroms)>0) warning(paste0("Some chromosomes, e.g. ",unmapped_chroms[1]," has no mapping,"))
	}

	mapping<-GRangesToMapping(
		ranges_to_map_to=where.to.map,
		chrom_suffix=chrom.suffix
	)
	#ranges	
	mapped<-unlist(liftOver(what.to.map,mapping$chain))
	#seqlength
	seqlengths(mapped)<-mapping$seqlengths
	return(mapped)
}


#' GRangesToMapping
#' 
#' GRangesToMapping creates an list that contain a \code{chain}, which is a \code{Chain} that descritbes a \code{liftOver} based on the intevals of the \code{ranges_to_map_to} parameter and a \code{seqlengths} that contain the lengths of the target (mapped) chromosomes. The mapping is from original chromosomes to the mapped chromosomes that are sticked intevals of the \code{ranges_to_map_to} per chromosome. 
#' 
#' @param ranges_to_map_to A \code{GRanges} file with non-overlapping intervals that will be converted to a chain file. Required.
#' @param chrom_suffix The suffix to be appended to all the sestination chromosome names in the mapping "default is "_mapped"
#' @param verbose Output updates while the function is running. Default FALSE
#' @return a list with two objects: \code{chain} is a \code{Chain} object for \code{liftOver} to map all the chromosomes according to GRanges; the \code{seqlegths} are the lenghthes of the mapped chromosomes
# this ia code by Veronica Busa and Alexander Favorov
GRangesToMapping<-function(ranges_to_map_to,
                                    chrom_suffix = "_mapped",
																		verbose=FALSE)
{
  #confirm GRanges doesn't have any overlapping intervals
  ranges_to_map_to<-reduce(ranges_to_map_to)
  #chromosomes to get lengths	
 	mapping<-list() 
	mapping$chain<-new('Chain')
  # make a chain for each chromosome in the genome
	mapping$seqlengths<-c()
	#Chain is an inner rtraclyaer class, here we try to explain it
	#it is list of ChainBlock -- one per chromosome in 'from' genome, the chromosome name is the index in the list,
	#names() give you all the chromosome names
	#each ChainBlick defines all the chains (in terms of liftover) that map from the chromosome
	#slot @ranges is a stack of all the ranges involved in the mapping
	#slot @offset is the shift in mapping for each interval, it is positive if the mapped space has higher position than the original
	#all the other slots are 'packed', I mean, thay could be given in vectors of the same length as the @ranges,
	#one item per range, and in this case the last slot, @length, is rep(1,lenght(@ranges)), 
	#and @score, @space and @revesed have the same length lenght(@ranges)
	#but actually, the equal values are joined, the length is how many times each value is repeated,
	#and @score, @space and @revesed have the same length lenght(@length)
	#each value - one per chain in terms of liftover
	#slot @score is chain score integer()
	#slot @space is where it maps (what chromosome in targed genome) character()
	#slot @reversed is boolean -- whether is goes reversed (complemntary chain) boolean()
	#slot @length is a vector od letghs of liftover chains (how many intervals are involde in each of them) integer()

	#in our simple case, each ChainBlock carries on chain, so all the last 4 slots are 1-element vectors 
  for(chr in as.character(ranges_to_map_to@seqinfo@seqnames)){
		#preparing ChainBlock, one per this chromosome
    if(verbose==TRUE){cat(paste("Chromosome", chr, "starting..."))}
    chr_ranges<-ranges(ranges_to_map_to %>% filter(seqnames==chr))
		len<-length(chr_ranges)
    if(len==0){cat(" empty \n");next;} # in case of chromosomes without data
		shift<-rep(as.integer(NA),len)
		#prepare the shifts
		shift[1]<-c(start(chr_ranges)[1]-1)
		#first shift -- now, the first mapped interval starts the mapped-to chromosome
		if (len>1) {
			for (i in 2:len) { #is is the index of interval inside the chromosome
					shift[i]<-shift[i-1]+start(chr_ranges)[i]-end(chr_ranges)[i-1]-1
					#accumutale left shift, the (start[i]-end[i-1]-1) left (e.g. negative) 
					#is what the i-th region acquired in addition to (i-1)-th 
			}
		}
		mapped_chr<-paste0(chr,chrom_suffix)
		mapping$chain[[chr]]<-new("ChainBlock",
			ranges=chr_ranges,
			offset=as.integer(shift),
			score=as.integer(c(42)),
			space=c(mapped_chr),
			reversed=c(FALSE),
			length=c(len)
		)
		mapping$seqlengths[mapped_chr]<-sum(width(chr_ranges))
	}
  return(mapping)
}

