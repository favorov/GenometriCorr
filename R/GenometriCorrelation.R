# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2014 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#
# MarkupsIndependence is the main function of the package

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
	if (inherits(namelist,"RangedData"))# if it is RangedData 
	{
		names(namelist)<-add.chr.prefix.to.names(names(namelist))
		return(namelist)
	}
	if (inherits(namelist,"GRanges"))
	{
		seqlevels(namelist)<-add.chr.prefix.to.names(seqlevels(namelist))
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

GenometriCorrelation <- function(
	query,reference,
	chromosomes.to.proceed=c(),
	chromosomes.to.include.in.awhole=c(),
	chromosomes.to.exclude.from.awhole=c(),
	add.chr.as.prefix=FALSE,
	awhole.only=FALSE,
	space=NA,
	map.to.half=TRUE,
	showProgressBar=TRUE,
	showTkProgressBar=FALSE,
	chromosomes.length=c(),
	suppress.evaluated.length.warning=FALSE,
	cut.all.over.length=FALSE,
	ecdf.area.permut.number=100,
	mean.distance.permut.number=100,
	jaccard.measure.permut.number=100,
	jaccard.permut.is.rearrangement=FALSE,
	awhole.space.name="awhole",
	keep.distributions=FALSE,
	representing.point.function=mitl,
	query.representing.point.function=representing.point.function,
	reference.representing.point.function=representing.point.function,
	supress.evaluated.length.warning
	)
{
	#query and reference are two data sets that are under analysis
	#each og them is IRanges or RangedData
	#list.of.spaces is list of spaces to work with
	#chromosomes.to.proceed is list of chromosomes to proceed
	#chromosomes.to.include.in.awhole is a list of chromosomes to include in 'whole genome', default is all chromosomes.to.proceed()
	#chromosomes.to.exclude.from.awhole excludes some chromosome from the list 
	#add.chr.as.prefix adds 'chr' prefix to all chromosome spaces if there is no 'Chr', 'chr', etc prefix at any of the names of the names array; all the prefixes like 'CHR', 'chr' etc are lowercased to 'chr'
	#space is the mane for one chromosome tested, it provides tha name if we compare IRanges,
	#or it takes one name from the list exactly like if give one name in chromosomes.to.proceed()
	#map.to.half is whether to calculate relative distances in [0,0.5] or in [0,1]
	#showProgressBar whether to show a progress indicator
	#showTkProgressBar whether to show a Tk progress indicator; work only if tcltk is loaded
	#chromosomes.length is an array of length of chromosomes, names are the names of chromosomes
	#cut.all.over.length if the length is given and there is any interval coord that is higher than the lenght, the default behavoiur (when cut.all.over.length is FALSE) is to show an error and stop. If it is TRUE, the interval will be just truncated.
	#suppress.evaluated.length.warning suppresses the warning that a chromosome lenght is eveluated rather then given. The evaluation is just the rightmost coord in all the intervals, it is used is the length is NA.
	#ecdf.area.permut.number is number of permutations for ecdf area method
	#mean.distance.permut.number is the same thing about the mean ref-to-query distance
	#awhole.space.name is the space name for this operation
	#keep.distributions if TRUE, the result list includes a the distributions that were used to obtain p-valies
	#representing.point.function is the function that calculates the representing point for any interval 
	#the default is mitl (calculates the middle),
	#the function is to look like mitl<-function(start,end,chromosome.length,space); the additional parameters that do not depend on chrom, etc are in ..
	#query.representing.point.function=representing.point.function for query
	#reference.representing.point.function=representing.point.function for reference
	#supress.evaluated.length.warning is a historical typo

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
		stop("The thing given as first (",nameQ,") range argument to\n  MarkupsIndependence is not an object!")
	if (!is.object(reference))  
		stop("The thing given as second (",nameR,") range argument to\n  MarkupsIndependence is not an object!")
	#they are both objects if we are here

	irQ<-irR<-rdQ<-rdR<-grQ<-grR<-FALSE #initialise

	if (!(irQ=inherits(query,"IRanges")) && !(rdQ=inherits(query,"RangedData")) && !(grQ=inherits(query,"GenomicRanges")))
		stop("The thing given as first (",nameQ,") range argument to MarkupsIndependence is\n  nor an IRanges nor a RangedData nor a GRanges object!")
	if (!(irR=inherits(reference,"IRanges")) && !(rdR=inherits(reference,"RangedData")) && !(grR=inherits(reference,"GenomicRanges")))
		stop("The thing given as second (",nameR,") range argument to MarkupsIndependence is\n  nor an IRanges nor a RangedData nor a GRanges object!")

	#they both are IRanges or RangedData if we are here

	#space is an old thing; currently, chromosomes.to.proceed is preferrable

	if (length(chromosomes.to.proceed)>0 && !is.na(space) && space != chromosomes.to.proceed[1])
		stop("Both chromosomes.to.proceed and space parameters are given to MarkupsIndependence.\nI do not know what to do")
	
	#if we are here, only one of chromosomes.to.proceed and space is given or chromosomes.to.proceed[1] == space
	if (is.na(space) && length(chromosomes.to.proceed))
		space<-chromosomes.to.proceed[1]
	
	if (!is.na(space) && length(chromosomes.to.proceed)==0)
		chromosomes.to.proceed[1]<-space

	if (is.na(space) && length(chromosomes.to.proceed)==0 &&
			length(chromosomes.length)==1 &&
			!is.null(names(chromosomes.length)))
			{
				space<-names(chromosomes.length)
				chromosomes.to.proceed<-c(space)
			}

	if ((length(chromosomes.length)==1) &&
			(
				is.null(names(chromosomes.length)) ||
				is.na(names(chromosomes.length)[1])
			))
			{
				names(chromosomes.length)<-c(space)
			}
	#now, chromosomes.to.proceed[1] and space are the same thing (possibly, space=="" and chromosomes.to.proceed is empty)
	

	if (add.chr.as.prefix)
	{
		space<-add.chr.prefix.to.names(space)
		chromosomes.to.proceed<-add.chr.prefix.to.names(chromosomes.to.proceed)
		chromosomes.to.exclude.from.awhole<-add.chr.prefix.to.names(chromosomes.to.exclude.from.awhole)
		chromosomes.to.include.in.awhole<-add.chr.prefix.to.names(chromosomes.to.include.in.awhole)
		names(chromosomes.length)<-add.chr.prefix.to.names(names(chromosomes.length))
		if(rdQ || grQ)
			query<-add.chr.prefix.to.names(query)
		if(rdR || grR)
			reference<-add.chr.prefix.to.names(reference)
	}

	#here, we define spacesA and spacesB	

	spacesA<-c()
	spacesB<-c()

	if (rdQ) spacesA<-names(query)
	if (grQ) spacesA<-seqlevels(query)
	if (rdR) spacesB<-names(reference)
	if (grR) spacesB<-seqlevels(reference)
	#list of space names as character vector

	if (irQ && !irR)
	{
		if(is.na(space))
		{
			if(length(spacesB)==1)
			{
				space<-spacesB[1]
				chromosomes.to.proceed[1]<-space
				query<-RangedData(space=c(space),ranges=query)
			}
			else
				stop("The first (",nameQ,") range argument is IRanges.\n",
				"  The second (",nameR,") range argument is RangedData with more than one chromosome (space).\n  The space argument is not defined.\n It's all lost!")
		}
		else
			query<-RangedData(space=c(space),ranges=query)
		spacesA<-spacesB<-c(space)
	}

	if (!irQ && irR)
	{
		if(is.na(space))
		{
			if(length(spacesA)==1)
			{
				space<-spacesA[1]
				chromosomes.to.proceed<-space
				reference<-RangedData(space=c(space),ranges=reference)
			}
			else
				stop("The first (",nameQ,") range argument is RangeData with more than one chromosome (space).\n",
				"  The second (",nameR,") range argument is IRanges.\n  The space argument is not defined.\n  It's all lost!")
		}
		else
			reference<-RangedData(space=c(space),ranges=reference)
		spacesA<-spacesB<-c(space)
	}
	
	if (irR && irQ)
	{
		if (space=="")	space<-"the_space"
		reference<-RangedData(space=c(space),ranges=reference)
		query<-RangedData(space=c(space),ranges=query)
		spacesA<-spacesB<-c(space)
	}

	# here, all the data is RangedData's

	list.of.spaces<-intersect(spacesA,spacesB)
	
	if (length(list.of.spaces)==0)
			stop("The query and reference have to have some common chromosome names, and they do not.\n",
			"  It's all lost!")

	#list.of.spaces now is what we can proceed

	#if chromosomes.to.proceed is empty, we proceed all
	if(length(chromosomes.to.proceed)==0) chromosomes.to.proceed<-list.of.spaces
	else chromosomes.to.proceed<-listtoString(chromosomes.to.proceed)

	#test prints
	#print(spacesA)
	#print(spacesB)
	#print(list.of.spaces)
	#print(chromosomes.to.proceed)

	list.of.spaces<-intersect(list.of.spaces,chromosomes.to.proceed)

	# now, list.of.spaces is what we will really proceed	
	if (length(list.of.spaces)==0)
		stop("There is no intersection between input set of common\n  chromosome names of query and reference and chromosomes.to.proceed.\n  It's all lost!")

	wrong_diff<-setdiff(chromosomes.to.proceed,list.of.spaces)
	if (length(wrong_diff)>0)
		warning("Some spaces are in chromosomes.to.proceed but they are not in the input data.\n They are:\n",paste(wrong_diff,collapse="\n"),"\n")

	#if the include list is empty, it is set to whole list
	if(length(chromosomes.to.include.in.awhole)==0) chromosomes.to.include.in.awhole<-list.of.spaces

	#we test whether all the chromosomes.to.include.in.awhole is in chromosomes.to.proceed	
	wrong_diff<-setdiff(chromosomes.to.include.in.awhole,list.of.spaces)
	if (length(wrong_diff)>0)
		stop("Some spaces  are in chromosomes.to.include.in.awhole but they are not in the process list.\n  It's all lost!\n",paste(wrong_diff,collapse="\n"), "\n")

	awhole.chromosomes<-intersect(chromosomes.to.include.in.awhole,list.of.spaces)

	wrong_diff<-setdiff(chromosomes.to.exclude.from.awhole,chromosomes.to.include.in.awhole)
	if (length(wrong_diff)>0)
		stop("Some spaces are in chromosomes.to.exclude.from.awhole but they are not in the process list.\n  It's all lost!:\n",paste(wrong_diff,collapse="\n"), "\n")

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

	if (grQ)
	{
		chromosomes.length.Q<-seqlengths(query)
		query<-as(query,'RangedData')
	}

	if (grR)
	{
		chromosomes.length.R<-seqlengths(reference)
		reference<-as(reference,'RangedData')
	}
	
	# we do not intereasted in lengths of those chromosomes that are
	# not in intersection of spacesA and spacesB;
	# if we have information in chromosomes.length, we do not care
	# about what was given in seqlengths of A and B

	for (name in list.of.spaces)
	{
		if (name %in% names(chromosomes.length) && 
				!is.na(chromosomes.length[name]) ) next
		lQ<-NA
		lR<-NA
		if (grQ) lQ<-chromosomes.length.Q[name]
		if (grR) lR<-chromosomes.length.R[name]

		if (!is.na(lQ) && !is.na(lR) && lQ!=lR)
		{
			warning("The lengths of chromosome", name ,"is different in two input GRanges,\n")
			next
		}
		if (!is.na(lQ)) chromosomes.length[name]<-lQ
		if (!is.na(lR)) chromosomes.length[name]<-lR
	}

	result<-.RangedDataGenometricsCorrelation(
			rd_query=query,
			rd_reference=reference,
			list.of.spaces=list.of.spaces,
			map.to.half=map.to.half,
			showProgressBar=showProgressBar,
			showTkProgressBar=showTkProgressBar,
			chromosomes.length=chromosomes.length,
			suppress.evaluated.length.warning=suppress.evaluated.length.warning,
			cut.all.over.length=cut.all.over.length,
			ecdf.area.permut.number=ecdf.area.permut.number,
			mean.distance.permut.number=mean.distance.permut.number,
			jaccard.measure.permut.number=jaccard.measure.permut.number,
			jaccard.permut.is.rearrangement=jaccard.permut.is.rearrangement,
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
	result@config$options$cut.all.over.length=cut.all.over.length
	result@config$tests=list()
	result@config$tests$ecdf.area.permut.number=ecdf.area.permut.number
	result@config$tests$mean.distance.permut.number=mean.distance.permut.number
	result@config$tests$jaccard.measure.permut.number=jaccard.measure.permut.number
	result@config$options$keep.distributions=keep.distributions
	if (awhole.space.name!="awhole") result@config$options$awhole.space.name=awhole.space.name

	#print(result@config)

	return(result)
}


.RangedDataGenometricsCorrelation<-function(
	rd_query,rd_reference,
	list.of.spaces,
	map.to.half=TRUE,
	showProgressBar=TRUE,
	showTkProgressBar=FALSE,
	chromosomes.length=c(),
	suppress.evaluated.length.warning=FALSE,
	cut.all.over.length=FALSE,
	ecdf.area.permut.number=100,
	mean.distance.permut.number=100,
	jaccard.measure.permut.number=100,
	jaccard.permut.is.rearrangement=FALSE,
	awhole.chromosomes=c(),
	awhole.space.name="awhole",
	awhole.only=FALSE,
	keep.distributions=FALSE,
	query.representing.point.function,
	reference.representing.point.function
	)
{

	#the thing actually calculates everything
	#it is to called from MarkupsIndependence
	#rd_query,rd_reference are two RangedDatas that are under analysis
	#list.of.spaces is list of spaces to work with
	#map.to.half is whether to calculate relative distances in [0,0.5] or in [0,1]
	#showProgressBar whether to show a progress indicator
	#showTkProgressBar whether to show a Tk progress indicator; work only if tcltk is loaded
	#chromosomes.length is an array of length of chromosomes, names are the names of chromosomes
	#ecdf.area.permut.number is number of permutations for ecdf area method
	#mean.distance.permut.number is the same thing about the mean ref-to-query distance
	#awhole.chromosomes is list of chromosomes to be included in whole-genome calculation
	#awhole.space.name is the space name for this operation
	
	if (length(awhole.chromosomes) > 0)
		do.awhole<-TRUE
	else 
		do.awhole<-FALSE
	
	# HJ -- tcltk now part of R; we can Depends it
	#if ((showTkProgressBar) && !require("tcltk",quietly=TRUE))
		#showTkProgressBar=FALSE

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
		if (! space %in% names(chromosomes.length))
			chromosomes.length[space]=NA
		if ( is.na(chromosomes.length[space]) ) #if it was absent, now it is NA as well as if it was NA
		{
			chromosomes.length[space]=chromosomes.length.eval(rd_query,rd_reference)
			if (! (suppress.evaluated.length.warning))
				warning(paste("Length for chromosome",space,"is evaluated rather than pre-given."))
		}
		else
		{
			if (sum(end( (ranges(rd_query))[[space]] ) > chromosomes.length[space]))
			{
				if (!cut.all.over.length)
				{
					coord=end( (ranges(rd_query))[[space]] )[(end( (ranges(rd_query))[[space]] ) > chromosomes.length[space])][[1]]
					stop('There is a query range (it ends at ',coord,') in chromosome ',space,' that spans out of the chromosome' )
				}

			}
			if (sum(end( (ranges(rd_reference))[[space]] ) > chromosomes.length[space]))
			{
				if (!cut.all.over.length)
				{
					coord=end( (ranges(rd_reference))[[space]] )[(end( (ranges(rd_reference))[[space]] ) > chromosomes.length[space])][[1]]
					stop('There is a query range (it ends at ',coord,') in chromosome ',space,' that spans out of the chromosome' )
				}
			}
		}
	}
	
	if (cut.all.over.length)
	{
		dA<-as.data.frame(rd_query)
		dB<-as.data.frame(rd_reference)
		#converted to dataframe
		used_names<-names(chromosomes.length)
		dA<-dA[as.character(dA[,'space']) %in% used_names,]
		dB<-dB[as.character(dB[,'space']) %in% used_names,]
		#removed all the spaces that ar not in names(chromosomes.length)
		dA<-dA[dA['start']<=as.vector(chromosomes.length[as.character(dA[,'space'])]),]
		dB<-dB[dB['start']<=as.vector(chromosomes.length[as.character(dB[,'space'])]),]
		#remove if start > length of the corresponging chromosome	
		#set end; for is easier :)
		for (row in 1:dim(dA)[1])
			if(dA[row,'end']>chromosomes.length[as.character(dA[row,'space'])])
				dA[row,'end']=chromosomes.length[as.character(dA[row,'space'])]
		for (row in 1:dim(dB)[1])
			if(dB[row,'end']>chromosomes.length[as.character(dB[row,'space'])])
				dB[row,'end']=chromosomes.length[as.character(dB[row,'space'])]
		rd_query<-RangedData(dA)
		rd_reference<-RangedData(dB)
	}

	good_space=rep(TRUE,length(list.of.spaces))

	for (space_no in 1:length(list.of.spaces)) #calculate everything nonpermutted for each chromosomes
	{
		space<-list.of.spaces[space_no]
		if (length(start(ranges(rd_query)[[space]]))==0)
			good_space[space_no]=FALSE
		if (length(start(ranges(rd_reference)[[space]]))==0)
			good_space[space_no]=FALSE
		if (!good_space[space_no])
			warning("The space: ",space, " is not populated, we omit it.")
	}

	list.of.spaces<-list.of.spaces[good_space];

	#filtered....

	#initialise it all
	result<-list()

	for (space in list.of.spaces)
	{
		result[[space]]<-list()
	}

	if (do.awhole)
	{
		result[[awhole.space.name]]<-list()

		result[[awhole.space.name]][['query.population']]<-0
		result[[awhole.space.name]][['reference.population']]<-0
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
			ranges=rd_query[space]$ranges,
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
			ranges=rd_reference[space]$ranges,
			representing.point.function=reference.representing.point.function,
			chromosome.length=chromosomes.length[space],
			space=space		)

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}

		result[[space]][['query.population']]<-length(qu)

		result[[space]][['reference.population']]<-length(ref)

	  result[[space]][['query.coverage']]<-sum(width(reduce(rd_query[space]$ranges)))

		result[[space]][['reference.coverage']]<-sum(width(reduce(rd_reference[space]$ranges)))
		
		result[[space]][['relative.distances.data']]<-
			query_to_ref_relative_distances(qu,ref,map.to.half,is_query_sorted=T,is_ref_sorted=T,chrom_length=chromosomes.length[space])

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}

		result[[space]][['relative.distances.ks.p.value']]<-
			ks.test(untie(result[[space]]$relative.distances.data),punif,min=0,max=rel.dist.top)$p.value

		if (showProgressBar) setTxtProgressBar(txt_pb, getTxtProgressBar(txt_pb)[1]+1)

		if (showTkProgressBar)
		{
			done <- getTkProgressBar(tk_pb)[1]+1;
			done_info <- sprintf('chromosome: %s ; %i of %i done',space,done,pb_capacity)
			setTkProgressBar(tk_pb,value=done,title=paste('GenometriCorrelation:',done_info),label=done_info)
		}
		
		result[[space]][['relative.distances.ecdf.deviation.area']]<-
			integrate(
					function(x){return(abs((ecdf(result[[space]]$relative.distances.data))(x)-x/rel.dist.top))},
					lower=0,upper=rel.dist.top,
					subdivisions=length(result[[space]]$relative.distances.data)*100,
					rel.tol=integr_rel_tol)$value
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

		result[[space]][['projection.test']]<-query.to.ref.projection.statistics(qu,rd_reference[space]$ranges,TRUE,chromosomes.length[space])	

		result[[space]][['query.reference.intersection']]<-sum(width(reduce(intersect(rd_query[space]$ranges,rd_reference[space]$ranges))))
		result[[space]][['query.reference.union']]<-sum(width(reduce(union(rd_query[space]$ranges,rd_reference[space]$ranges))))

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

		proj.p.value<-
			pbinom(
				result[[space]][['projection.test']][['query.hits']],
				result[[space]][['query.population']],
				result[[space]][['projection.test']][['reference.coverage']]/result[[space]][['projection.test']][['space.length']]
			)

		if (proj.p.value<.5) # one-sided test
		{
			result[[space]][['projection.test.p.value']]<-proj.p.value
			result[[space]][['projection.test.lower.tail']] <- TRUE
		}
		else
		{
			result[[space]][['projection.test.p.value']]<- 1.-proj.p.value
			result[[space]][['projection.test.lower.tail']] <- FALSE
		}

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
		the.names<- names(result[[awhole.space.name]])
		if ('relative.distances.data' %in% the.names)
		{
			space<-awhole.space.name	
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
			space<-awhole.space.name	
			proj.p.value<-
				pbinom(
					result[[space]][['projection.test']][['query.hits']],
					result[[space]][['query.population']],
					result[[space]][['projection.test']][['reference.coverage']]/result[[space]][['projection.test']][['space.length']]
				)
			if (proj.p.value<.5) # one-sided test
			{
				result[[space]][['projection.test.p.value']]<-proj.p.value
				result[[space]][['projection.test.lower.tail']] <- TRUE
			}
			else
			{
				result[[space]][['projection.test.p.value']]<- 1.-proj.p.value
				result[[space]][['projection.test.lower.tail']] <- FALSE
			}
		}

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
				sample_size<-length(rd_query[space]$ranges)
				sample<-runif(sample_size,0,rel.dist.top) #sample
				dev_area<-
					integrate(
							function(x){return(abs((ecdf(sample)(x))-x/rel.dist.top))},
							lower=0,upper=rel.dist.top,
							subdivisions=sample_size*100,
							rel.tol=integr_rel_tol)$value

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
				sample_size<-length(rd_query[space]$ranges)
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
					permuted.query<-.permuteIRanges(rd_query[space]$ranges,chr_length)
				else
					permuted.query<-.rearrangeIRanges(rd_query[space]$ranges,chr_length)

				space.intersection<-sum(width(reduce(intersect(permuted.query,rd_reference[space]$ranges))))
				space.union<-sum(width(reduce(union(permuted.query,rd_reference[space]$ranges))))
				
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
			p.value<-
				1-(ecdf(result[[space]][['relative.distances.ecdf.deviation.area.null.list']]))(result[[space]][['relative.distances.ecdf.deviation.area']]) #we treat only right side
			if (p.value==0) 
				p.value = paste("<",toString(1/ecdf.area.permut.number),sep='')
			result[[space]][['relative.distances.ecdf.deviation.area.p.value']]<-p.value
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
			p.value<-
				(ecdf(result[[space]][['scaled.absolute.min.distance.sum.null.list']]))(result[[space]][['scaled.absolute.min.distance.sum']]) #we treat only right side
			
			if(p.value<=0.5)
				lower.tail<-TRUE
			else
			{
				lower.tail <- FALSE
				p.value<-1-p.value
			}

			if (p.value==0) 
				p.value = paste("<",toString(1/mean.distance.permut.number),sep='')

			result[[space]][['scaled.absolute.min.distance.sum.p.value']]<-p.value
			result[[space]][['scaled.absolute.min.distance.sum.lower.tail']] <- lower.tail 
		}
	}
	#jaccard distance sum p-values
	if (jaccard.measure.permut.number>0)
	{
		for (space in c(list.of.spaces,awhole.space.name))
		{
			p.value<-
				(ecdf(result[[space]][['jaccard.measure.null.list']]))(result[[space]][['jaccard.measure']]) #we treat only right side
			
			if(p.value<=0.5)
				lower.tail<-TRUE
			else
			{
				lower.tail <- FALSE
				p.value<-1-p.value
			}

			if (p.value==0) 
				p.value = paste("<",toString(1/jaccard.measure.permut.number),sep='')

			result[[space]][['jaccard.measure.p.value']]<-p.value
			result[[space]][['jaccard.measure.lower.tail']] <- lower.tail 
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


sorted.representing.points<-function(ranges,representing.point.function,chromosome.length,space)
{
	#ranges is IRanges
	mids<-representing.point.function(start(ranges),end(ranges),chromosome.length,space)
	return(mids[order(mids)])
}

#input: Irange
#output: vector of represetating points


chromosomes.length.eval<-function(rd_query,rd_reference)
{
	#rd_query and rd_reference are IRange
	min_coord<-min(start(rd_query),start(rd_reference),end(rd_query),end(rd_reference),0)
	#we prefer this way rather than obvious 0
	#to equal the error if everything is located in say right and left telomeres
	max_coord<-max(start(rd_query),start(rd_reference),end(rd_query),end(rd_reference))
	return(max_coord-min_coord)
}


query_to_ref_relative_distances<-function(query,ref,map.to.half,is_query_sorted=F,is_ref_sorted=F,chrom_length=NA)
#calculate distances for a pair of positions (points)  vectors: query, ref
{

	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]
	rel_data<-list()

	#if chrom_length==0, we will loose all probes that fall out ref range, 
	#otherwise we provide the chrom length and so circle the reference
	refindex<-1
	firstref<-ref[1]
	lastref<-ref[length(ref)]
	rel_data<-c()
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

	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]
	rel_data<-list()

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

	if (!is_query_sorted) 
		query<-query[order(query)]
	if (!is_ref_sorted) 
		ref<-ref[order(ref)]
	rel_data<-list()

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

	if (!is_query_sorted) 
		query<-query[order(query)]

	#we reduce it, whatever!
	ref<-reduce(ref) #remove duplication in coverage

	ref<-ref[order(ref)]
	

	#if chrom_length==0, we will loose all probes that fall out ref range, 
	#otherwise we provide the chrom length and so circle the reference
	currentref<-1
#	nextcurrentref<-1
	projection_data<-c()
	projection_data[['space.length']]<-chrom_length
	projection_data[['reference.coverage']]<-sum(width(reduce(ref)))
	
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

