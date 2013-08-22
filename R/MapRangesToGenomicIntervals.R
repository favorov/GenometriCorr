# GenometriCorrelation project evaluating two markups genomwide independence. 
# (c) 2010-2011 Alexander Favorov, Leslie Cope, Yulia Medvedeva, 
#              Loris Mularoni, Vsevolod Makeev, Sarah Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors
# $Id: MapRangesToGenomicIntervals.R 1717 2012-04-18 19:29:15Z favorov $



MapRangesToGenomicIntervals<-function(
	where.to.map, what.to.map,
	chromosomes.to.proceed=NA,
	unmapped.chromosome.warning=TRUE,
	unmapped.range.warning=FALSE
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
	ranges<-IRanges()
	seqnames<-c()
	seqleninfo<-c()
	mapped_start=c()
	mapped_end=c()
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
		if (! isNormal(IRanges(start=start(whereranges),end=end(whereranges)-1)))
		{
			warning(paste("The chromosome",chr,"is not normalised in where.to.map; skipped."))
			next;
		}
		whatranges<-sort(ranges(what.to.map)[[chr]])
		#we decide what to do with the chromosome lengths
		ind_where<-0
		ind_what<-0
		#print('What')
		#print(whatranges)
		#print('Where')
		#print(whereranges)
		while (ind_where < length(whereranges))
		{
			ind_where<-ind_where+1
			mapped_range_name=paste(chr,start(whereranges)[ind_where],sep='_');
			#cat(chr,mapped_range_name,ind_where,"\n")
			this_pseudocromosome_population<-0
			while((ind_what < length(whatranges)) && (start(whatranges)[ind_what+1] < end(whereranges)[ind_where]))
			{
				ind_what<-ind_what+1
				if ( end(whatranges)[ind_what] < start(whereranges)[ind_where] ) next; 
				#moving 'what' cursor up to get current 'where' interval
				cur_start<-max(1,start(whatranges)[ind_what]-start(whereranges)[ind_where]+1)
				cur_end<-min(end(whatranges)[ind_what]-start(whereranges)[ind_where]+1,width(whereranges)[ind_where])
				mapped_start<-c(mapped_start,cur_start)
				mapped_end<-c(mapped_end,cur_end)
				this_pseudocromosome_population<-this_pseudocromosome_population+1
				#cat('chr:',chr,'\n')
				#cat('wha: ',ind_what,':',start(whatranges)[ind_what],end(whatranges)[ind_what],"\n")
				#cat('whe: ',ind_where,':',start(whereranges)[ind_where],end(whereranges)[ind_where],"\n")
				#cat('res: ',length(mapped_end),':',cur_start,cur_end,"\n")
			}
			seqleninfo<-c(seqleninfo,width(whereranges)[ind_where])
			names(seqleninfo)[[length(seqleninfo)]]<-mapped_range_name
			seqnames<-c(seqnames,rep(mapped_range_name,this_pseudocromosome_population))
		}
	}
	if(unmapped.range.warning && ( length(space(what.to.map)) > length(mapped_start) ))
		warning("Some ranges remained unmapped.\n")
	if (length(seqnames)==0)
	{
		warning("Empty mapping result.\n");
		return (GRanges())
	}
	return
	(
		GRanges(
			seqnames=seqnames,
			ranges=IRanges(start=mapped_start,end=mapped_end),
			seqlengths=seqleninfo)
	)
}

#isNormal(ranges(cpgis)[as.vector(unique(space(cpgis)))])
