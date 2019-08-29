# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2018 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors



#'@export
MapRangesToGenomicIntervals<-function(
	where.to.map, what.to.map,
	chromosomes.to.proceed=NA,
	unmapped.chromosome.warning=TRUE,
	unmapped.range.warning=FALSE,
	nonnormalised.mapping.warning=TRUE
)
#subgenome and rd are suppose to be RangedData or GRanges
#in the second case, we test the lenthg equivalence
{
	is.where.gr<-inherits(where.to.map,"GRanges")
	is.where.rd<-inherits(where.to.map,"RangedData")
	if (!is.where.rd && !is.where.gr)
		stop("where.to.map is not RangedData and it is not GRanges. It's all lost!\n")	
	chromosome.names.where<-NA
	if (is.where.gr)
	{
		chromosome.names.where<-as.character(unique(seqnames(where.to.map)))
		where.to.map<-as(where.to.map,'RangedData')
	}
	else #is.where.rd
	{
		chromosome.names.where<-as.character(unique(space(where.to.map)))
	}

	is.what.gr<-inherits(what.to.map,"GRanges")
	is.what.rd<-inherits(what.to.map,"RangedData")
	if (!is.what.rd && !is.what.gr)
		stop("what.to.map is not RangedData and it is not GRanges. It's all lost!\n")	
	chromosome.names.what<-NA
	if (is.what.gr)
	{
		chromosome.names.what<-as.character(unique(seqnames(what.to.map)))
		what.to.map<-as(what.to.map,'RangedData')
	}
	else #is.what.rd
	{
		chromosome.names.what<-as.character(unique(space(what.to.map)))
	}
	seqnames<-c()
	seqleninfo<-c()
	start=c()
	end=c()
	for (chr in chromosome.names.what)
	{
		if (!is.na(chromosomes.to.proceed) && ! chr %in% chromosomes.to.proceed) next;
		if (! chr %in% chromosome.names.where )
		{
			if (unmapped.chromosome.warning) 
				warning(paste("The chromosome",chr,"has no space to be mapped to."))
			next;
		}
		whereranges<-sort(ranges(where.to.map)[[chr]])

		if (! isNormal(whereranges))
		{
			if (nonnormalised.mapping.warning) 
				warning(paste("The chromosome",chr,"is not normalised in where.to.map; I do it for you."))
			whereranges <- asNormalIRanges(whereranges)
		}
	
		whatranges<-sort(ranges(what.to.map)[[chr]])
		
		mapping.mat<-as.matrix(findOverlaps(whatranges,whereranges))
		
		map.result<-apply(mapping.mat,1,
			function(hit){
				ind.what<-hit[1]
				ind.where<-hit[2]
				start<-max(1,start(whatranges)[ind.what]-start(whereranges)[ind.where]+1)
				end<-min(end(whatranges)[ind.what]-start(whereranges)[ind.where]+1,width(whereranges)[ind.where])
				seqnames<-paste0(chr,':',start(whereranges)[ind.where],'-',end(whereranges)[ind.where])
				c(seqnames=seqnames,end=end,start=start)
			}
		)
	
		seqnames<-c(seqnames,map.result['seqnames',])
		start<-c(start,as.integer(map.result['start',]))
		end<-c(end,as.integer(map.result['end',]))
	}
	if(unmapped.range.warning && ( length(space(what.to.map)) > length(start) ))
		warning("Some ranges remained unmapped.\n")
	if (length(seqnames)==0)
	{
		warning("Empty mapping result.\n");
		return (GRanges())
	}
	#prepare seqlength info
	seqlenames<-unique(seqnames)
	seqlengths<-sapply(as.vector(strsplit(seqlenames,':')), #list, we need
		function(splt){
			startend<-as.integer(strsplit(splt[2],'-')[[1]])
			startend[2]-startend[1]+1	
		}
	)
	names(seqlengths)<-seqlenames
	return
	(
		GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=start,end=end),
			seqlengths=seqlengths)
	)
}

