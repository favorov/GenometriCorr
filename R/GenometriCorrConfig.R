# GenometriCorrelation project evaluating two genometric annotations 
#	correlation. Also provides some service IRanges-related procedures. 
# (c) 2010-2011 Alexander Favorov, Leslie Cope, Yulia Medvedeva, 
#              Loris Mularoni, Vsevolod Makeev, Sarah Wheelan.
#
# configuration - works with ini file format  
# $Id: GenometriCorrConfig.R 1893 2013-07-08 08:10:12Z favorov $

#.Parse.INI part of code is by: 
#Earl F. Glynn
#Scientific Programmer
#Stowers Institute for Medical Research
#https://stat.ethz.ch/pipermail/r-help/2007-June/134115.html
# Prototype of how to read INI files to process olfactometer data
# efg, 13 June 2007
# Thanks to Gabor Grothendieck for helpful suggestions in the R-Help
# mailing list on how to parse the INI file.

#if (!require('methods')) stop('GenometriCorrConfig requires methods package!\n')

.Parse.INI <- function(INI.filename)
{
  connection <- file(INI.filename)
  Lines  <- readLines(connection)
  close(connection)

  Lines <- chartr("[]", "==", Lines)  # change section headers

  connection <- textConnection(Lines)
  d <- read.table(connection, as.is = TRUE, sep = "=", fill = TRUE)
  close(connection)

	
	L <- d$V1 == ""                    # location of section breaks
  
	#d <- transform(d, V3=V2[which(L)[cumsum(L)]])[1:3] 

	d$V3=d$V2[which(L)[cumsum(L)]] #the same; R CMD check likes it more

	#L is 1 for sect header; cumsum(L) is # if section we are in; 
	#which(L)[cumsum(L)] shows where to see the header
 	
	#d <- subset(d, V1!="")

	d<-d[!L,] #the same; R CMD check likes it more


	d$V2[d$V2=='ON']='TRUE'
	d$V2[d$V2=='YES']='TRUE'
	d$V2[d$V2=='OFF']='FALSE'
	d$V2[d$V2=='NO']='FALSE'

  ToParse  <- paste("INI.list$", d$V3, "$",  d$V1, " <- '",
                    d$V2, "'", sep="")

  INI.list <- list()
  eval(parse(text=ToParse))

	#converting string values to values...

	if (!is.null(INI.list$data$do.mapping)) INI.list$data$do.mapping=eval(parse(text=INI.list$data$do.mapping))

	section='options'
	for (name in names(INI.list[[section]]))
		INI.list[[section]][[name]]=eval(parse(text=INI.list[[section]][[name]]))

	section='tests'
	for (name in names(INI.list[[section]]))
		INI.list[[section]][[name]]=eval(parse(text=INI.list[[section]][[name]]))

	section='chromosomes.length'
	for (name in names(INI.list[[section]]))
		INI.list[[section]][[name]]=eval(parse(text=INI.list[[section]][[name]]))

	if (is.null(INI.list$data) && !is.null(INI.list$files)) 
	#it is an old-fashioned file for <=1.0.15
	#it has [files] intead of [data]
		names(INI.list)[names(INI.list)=='files']='data'

	return(INI.list)
}


setClass('GenometriCorrConfig',contains='namedList',representation(src='character'))

setMethod('initialize', 'GenometriCorrConfig', function(.Object, src="")
	{
		if (src!="")
		{
			.Object@src=src
			.Object@.Data<-.Parse.INI(src)
			#for (name in names(ini))
			#	.Object@data[name]<-ini[name]
		}
		else #default
		{
			.Object@src=''
			.Object@.Data<-list()
		}
		.Object
	})

setMethod('show','GenometriCorrConfig',function(object)
	{
		for (name in names(object))
		{
			cat('[',name,']\n',sep='')
			section<-(object)[[name]]
			for( tag in names(section))
				if(!is.null(section[[tag]][1]) && (tagged<-section[[tag]][1])!="")
					cat(tag,'=',tagged,'\n',sep='')
				else
					cat(tag,'\n',sep='')
			cat('\n')				
		}
	})


setGeneric('run.config',
	function(
		conf,
		query=NA,
		reference=NA,
		query.format=NA,
		reference.format=NA,
		mapping=NA,
		mapping.format=NA,
		do.mapping=NA
		) standardGeneric('run.config'))

setMethod('run.config', signature(conf='GenometriCorrConfig'),
	function(
		conf,
		query=NA,
		reference=NA,
		query.format=NA,
		reference.format=NA,
		mapping=NA,
		mapping.format=NA,
		do.mapping=NA)
	{
		#[data] :
		#if the query.format parameter the value 'R.variable.name' 
		#we take the object with the name 'query' from the caller environment
		#any other query.format value (including the default NA) suggests a 
		#text file of the format to be read if query is character
		#if query (parameter, it cannot be in config file)
		#is a IRanges/RangedData/GRanges object, the run.config
		#uses the object and it does not write the
		#result@config lines for query and query.format
		#the name of the file is given by query variable
		#the same for reference and for mapping

		query.as.is<-reference.as.is<-mapping.as.is<-FALSE

		#object cannot be NA
		if (is.object(query) && 
				(inherits(query,"RangedData") || inherits(query,'IRanges') || inherits (query,'GRanges'))
			) 
		{
			#if query is an object, we do not go on with the $data stuff
			query.as.is<-TRUE
			conf$data[['query']]<-NULL
			conf$data[['query.format']]<-NULL
		}
		else if (!is.na(query))
		{
			#it is not an object, we want to test whether it is a string
			if (!is.na(query) && !is.character(query) )
				stop("The query parameter is not character and not an IRanges or RangedData or GRanges.\n")
			if (!is.na(query.format) && !is.character(query.format) )
				stop("The query.format parameter is not character.\n")
			conf$data$query<-query
			if (! is.na(query.format) ) conf$data$query.format<-query.format
		}
		else
		{
			if (is.null(conf$data[['query']])) # we need [[]] for exact match
				stop("No query parameter!\n")
		}

		if ((! is.null(conf$data$query.format) && conf$data$query.format=='R.variable.name'))
		# query is name of a variable.... Claculating....
		{
			#test whether it is a real var or stop
			query.name<-conf$data$query
			query<-NA
			query<-try(eval(parse(text=query.name)),silent=TRUE)
			if(!is.object(query))
			{
				if (is.na(query))
					stop(paste("The (query) variable with name:",query.name,"does not exist.\n"))
				else
					stop(paste("The (query) variable with name:",query.name,"is not an object.\n"))
			}		
			else if ( !( inherits(query,"RangedData") || inherits(query,'IRanges') || inherits(query,'GRanges') ) )
				stop(paste("The (query) variable with name:",query.name,"is not an IRanges or RangedData or GRanges.\n"))
			query.as.is<-TRUE	
		}


		if(!query.as.is)
		{
			if(!require("rtracklayer",quietly=TRUE))
			{
				source("http://bioconductor.org/biocLite.R")
				biocLite("rtracklayer")
				require("rtracklayer")
			}

			if (is.null(conf$data$query) || conf$data$query == '') stop ('Empty query filename\n')
			if (! file.exists(conf$data$query)) stop(paste('Cannot read from query filename: ',conf$data$query))

			if (is.null(conf$data$query.format)) 
				query<-rtracklayer::import(conf$data$query)
			else if (conf$data$query.format=='bed.like.with.header')
				query<-readTableToIRanges(conf$data$query,header=TRUE)
			else query<-rtracklayer::import(conf$data$query,format=conf$data$query.format)
		}

		#object cannot be NA
		if (is.object(reference) && 
				(inherits(reference,"RangedData") || inherits(reference,'IRanges') || inherits (reference,'GRanges'))
			) 
		{
			#if reference is an object, we do not go on with the $data stuff
			reference.as.is<-TRUE
			conf$data[['reference']]<-NULL
			conf$data[['reference.format']]<-NULL
		}
		else if (!is.na(reference))
		{
			#it is not an object, we want to test whether it is a string
			if (!is.na(reference) && !is.character(reference) )
				stop("The reference parameter is not character and not an IRanges or RangedData or GRanges.\n")
			if (!is.na(reference.format) && !is.character(reference.format) )
				stop("The reference.format parameter is not character.\n")
			conf$data$reference<-reference
			if (! is.na(reference.format) ) conf$data$reference.format<-reference.format
		}
		else
		{
			if (is.null(conf$data[['reference']])) # we need [[]] for exact match
				stop("No reference parameter!\n")
		}

		if ((! is.null(conf$data$reference.format) && conf$data$reference.format=='R.variable.name'))
		# reference is name of a variable.... Claculating....
		{
			#test whether it is a real var or stop
			reference.name<-conf$data$reference
			reference<-NA
			reference<-try(eval(parse(text=reference.name)),silent=TRUE)
			if(!is.object(reference))
			{
				if (is.na(reference))
					stop(paste("The (reference) variable with name:",reference.name,"does not exist.\n"))
				else
					stop(paste("The (reference) variable with name:",reference.name,"is not an object.\n"))
			}		
			else if ( !( inherits(reference,"RangedData") || inherits(reference,'IRanges') || inherits(reference,'GRanges') ) )
				stop(paste("The (reference) variable with name:",reference.name,"is not an IRanges or RangedData or GRanges.\n"))
			reference.as.is<-TRUE	
		}


		if(!reference.as.is)
		{
			if(!require("rtracklayer",quietly=TRUE))
			{
				source("http://bioconductor.org/biocLite.R")
				biocLite("rtracklayer")
				require("rtracklayer")
			}

			if (is.null(conf$data$reference) || conf$data$reference == '') stop ('Empty reference filename\n')
			if (! file.exists(conf$data$reference)) stop(paste('Cannot read from reference filename: ',conf$data$reference))


			if (is.null(conf$data$reference.format)) reference<-rtracklayer::import(conf$data$reference)
			else if (conf$data$reference.format=='bed.like.with.header')
				reference<-readTableToIRanges(conf$data$reference,header=TRUE)
			else reference<-rtracklayer::import(conf$data$reference,format=conf$data$reference.format)
		}

		#default is false; paramter overrides file; if default value is used, we do not 
		# show it in result@config
		if (!is.na(do.mapping)) conf$data$do.mapping<-do.mapping
		if (!is.null(conf$data$do.mapping)) do.mapping<-conf$data$do.mapping
		else do.mapping<-FALSE



		{#old
			#if (! is.na(mapping) ) conf$data$mapping<-mapping
			#if (! is.na(mapping.format) ) conf$data$mapping.format<-mapping.format
			#if (! is.na(mapping) && ( is.na(mapping.format) || mapping.format!='none')) conf$data$do.mapping<-TRUE
		}

		chr.len<-c()

		if (!is.null(conf$data$do.mapping) && conf$data$do.mapping)
		{

			#object cannot be NA
			if (is.object(mapping) && 
					(inherits(mapping,"RangedData") || inherits(mapping,'IRanges') || inherits (mapping,'GRanges'))
				) 
			{
				#if mapping is an object, we do not go on with the $data stuff
				mapping.as.is<-TRUE
				conf$data[['mapping']]<-NULL
				conf$data[['mapping.format']]<-NULL
			}
			else if (!is.na(mapping))
			{
				#it is not an object, we want to test whether it is a string
				if (!is.na(mapping) && !is.character(mapping) )
					stop("The mapping parameter is not character and not an IRanges or RangedData or GRanges.\n")
				if (!is.na(mapping.format) && !is.character(mapping.format) )
					stop("The mapping.format parameter is not character.\n")
				conf$data$mapping<-mapping
				if (! is.na(mapping.format) ) conf$data$mapping.format<-mapping.format
			}
			else
			{
				if (is.null(conf$data[['mapping']])) # we need [[]] for exact match
					stop("No mapping parameter!\n")
			}

			if ((! is.null(conf$data$mapping.format) && conf$data$mapping.format=='R.variable.name'))
			# mapping is name of a variable.... Claculating....
			{
				#test whether it is a real var or stop
				mapping.name<-conf$data$mapping
				mapping<-NA
				mapping<-try(eval(parse(text=mapping.name)),silent=TRUE)
				if(!is.object(mapping))
				{
					if (is.na(mapping))
						stop(paste("The (mapping) variable with name:",mapping.name,"does not exist.\n"))
					else
						stop(paste("The (mapping) variable with name:",mapping.name,"is not an object.\n"))
				}		
				else if ( !( inherits(mapping,"RangedData") || inherits(mapping,'IRanges') || inherits(mapping,'GRanges') ) )
					stop(paste("The (mapping) variable with name:",mapping.name,"is not an IRanges or RangedData or GRanges.\n"))
				mapping.as.is<-TRUE	
			}


			if(!mapping.as.is)
			{
				if(!require("rtracklayer",quietly=TRUE))
				{
					source("http://bioconductor.org/biocLite.R")
					biocLite("rtracklayer")
					require("rtracklayer")
				}

				if (is.null(conf$data$mapping) || conf$data$mapping == '') stop ('Empty mapping filename\n')
				if (! file.exists(conf$data$mapping)) stop(paste('Cannot read from mapping filename: ',conf$data$mapping))


				if (is.null(conf$data$mapping.format)) mapping<-rtracklayer::import(conf$data$mapping)
				else if (conf$data$mapping.format=='bed.like.with.header')
					mapping<-readTableToIRanges(conf$data$mapping,header=TRUE)
				else mapping<-rtracklayer::import(conf$data$mapping,format=conf$data$mapping.format)
			}

			chrom<-NA #default value means nothing
			if (!is.null(conf[['chromosomes']])) # not to mix with chromosomes.length
				chrom<-names(conf[['chromosomes']])

			query<-
				MapRangesToGenomicIntervals(
					where.to.map=mapping,
					what.to.map=query,
					unmapped.chromosome.warning=FALSE,
					chromosomes.to.proceed=chrom)
			reference<-
				MapRangesToGenomicIntervals(
					where.to.map=mapping,
					what.to.map=reference,
					unmapped.chromosome.warning=FALSE,
					chromosomes.to.proceed=chrom)
		}
		else
		{
			for (name in names(conf$chromosomes.length))
				chr.len[name]<-as.integer(conf$chromosomes.length[[name]][1])
		}# on mapping, no chromosomes.length is necessary

		if(!is.null(conf$tests$random.seed)) set.seed(conf$tests$random.seed)

		todo<-'result<-GenometriCorrelation(query=query,reference=reference'
		
		if ((is.null(conf$data$do.mapping) || !conf$data$do.mapping) &&	!is.null(conf[['chromosomes']]))
			todo<-paste(todo,',chromosomes.to.proceed=names(conf[[\'chromosomes\']])',sep='')
		#chromosomes list makes sense here only if mapping is off

		if (!is.null(conf$options$add.chr.as.prefix))
			todo<-paste(todo,',add.chr.as.prefix=',conf$options$add.chr.as.prefix,sep='')

		if (!is.null(conf$options$awhole.only))
			todo<-paste(todo,',awhole.only=',conf$options$awhole.only,sep='')
		
		if (!is.null(conf$options$showProgressBar))
			todo<-paste(todo,',showProgressBar=',conf$options$showProgressBar,sep='')

		if (!is.null(conf$options$showTkProgressBar))
			todo<-paste(todo,',showTkProgressBar=',conf$options$showTkProgressBar,sep='')

		if (length(chr.len) != 0) 
			todo<-paste(todo,',chromosomes.length=chr.len',sep='')

		if (!is.null(conf$options$supress.evaluated.length.warning))
			todo<-paste(todo,',supress.evaluated.length.warning=',conf$options$supress.evaluated.length.warning,sep='')

		if (!is.null(conf$options$cut.all.over.length))
			todo<-paste(todo,',cut.all.over.length=',conf$options$cut.all.over.length,sep='')

		if (!is.null(conf$tests$jaccard.measure.permut.number))
			todo<-paste(todo,',jaccard.measure.permut.number=',conf$tests$jaccard.measure.permut.number,sep='')

		if (!is.null(conf$tests$mean.distance.permut.number))
			todo<-paste(todo,',mean.distance.permut.number=',conf$tests$mean.distance.permut.number,sep='')

		if (!is.null(conf$tests$ecdf.area.permut.number))
			todo<-paste(todo,',ecdf.area.permut.number=',conf$tests$ecdf.area.permut.number,sep='')

		if (!is.null(conf$options$keep.distributions))
			todo<-paste(todo,',keep.distributions=',conf$options$keep.distributions,sep='')

		todo<-paste(todo,')',sep='')

		eval(parse(text=todo)) #result<-GenometriCorrelation(parameters)

		result@config<-conf # GenometriCorrelation write conf with default values, here we have one already, and we show it in result@config

		return(result)
	})


#setMethod('run.config',signature(conf='GenometriCorrConfig'),function(conf)
#	{
#	})
#run.config<-function(conf)
