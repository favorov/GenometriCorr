# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2016 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
# readTableToIRanges.R reads a variety of GFF- or BED-like  formats to an IRanges or a RangedData object


## file is the name of file to read
##
## header indicates whether the first noncomment row is interpreted as a header row
## if header is TRUE, each column is named by the its label in the header raw
##
## space is the number (1-based) or the name in the header of the column containing the 
##			chromosome name in the RangedData notation, it is 'space'. 
##      If space is NA or is not given, the data is read to IRanges object
##			unless any additional features (see below) are claimed to be read in the object
##
## start is the number or the name in the header of column with the interval start
##
## end is the number or the name in the header of column with the interval end
##
## width is the number or the name in the header of column with the interval width
##
## sep is the separator used in the file. "" means a whitespace-delitited table
##
## skip is number of omitted lines in the head of the file
##
## comment.char is a prefix of lines to be ignored
##
## autofeatures works if header is TRUE. It adds all the columns that are not interpreted as space, end, start or width as the
##      features of the intervals to the resulting data 
##
## ... contains parameters in form: feature_name=col_number or feature_name=column name, name is not one of the decribed parameter list.
##			Each of these pairs will add the feature feature_name that is read from the file's column to the resulting data
##
##      If autofeatures is TRUE or any ... pair is given the return is RangedData; 
##      if space is NA or not given, all the inetrvals are read in the same space;
##      the name of the space is the file name (see first line)
##



readTableToIRanges <- function(file=NA, space=NA, start=NA, end=NA, width=NA, sep="", skip=0, comment.char="#",header=F,autofeatures=F,...)
{
	max_chrom<-10000

	if (is.na('file'))
		stop("The file name is missing.");

	if (!is.character(file))
		stop("The file name is not a string.");

	if (autofeatures && !header)
		warning('The autofeatures switch is ignored because the header is off.')

	#test header to be a boolean
	
	#read the table
	q<-read.table(file=file,sep=sep,skip=skip,comment.char=comment.char,header=header)  #,stringsAsFactors=FALSE

	num_col=dim(q)[2]

	if (num_col <= 2) stop(paste("The file has less that 2 columns, maybe the separator sep=\"",sep,"\" is wrong?",sep=""))

	
	#print (c("le=",length(q),"           ",di=dim(q)))
	if (header) 
	{
		the_names<-names(q)
		columns_accounted=rep(F,length(the_names))
	}
	else
	{
		columns_accounted<-rep(FALSE,num_col)
	}
	
	space_names<-c('Chromosome','Chr','Space')

	if (is.na(space))
	{
		if (header)
		{
			for (name in space_names)
			{
				li<-grep(name,the_names,ignore.case=T)
				if (length(li)<1) next
				space<-li[1];
				break;
			}
		}
	}
	else
	{
		if(header)
		{
			li<-grep(space,the_names,ignore.case=T)
			if (length(li)>=1)
				space<-li[1];
		}
		if (is.na(as.integer(space)))
			stop("The space column is not an integer number nor a header name.")
		else if (space < 1)
			stop("The space column number is not a natural number.");
	}
	
	if (!is.na(space)) columns_accounted[space]=T

	start_names<-c('Start','Begin','First')

	
	if (is.na(start))
	{
		if (header)
		{
			for (name in start_names)
			{
				li<-grep(name,the_names,ignore.case=T)
				if (length(li)<1) next
				start<-li[1];
				break;
			}
		}
		if (is.na(start)) #did not help
			stop("The start column number is not given.")
	}
	else
	{	
		if(header)
		{
			li<-grep(start,the_names,ignore.case=T)
			if (length(li)>=1)
				start<-li[1];
		}
		if (is.na(start<-as.integer(start)))
			stop("The start column number is not an integer number nor a header name.")
		else if (start < 1)
			stop("The start column number is not a natural number.")
	}

	if (!is.na(start)) columns_accounted[start]=T

	end_names=c('End','Last','Stop')

	if (is.na(end))
	{
		if (header)
		{
			for (name in end_names)
			{
				li<-grep(name,the_names,ignore.case=T)
				if (length(li)<1) next
				end<-li[1];
				break;
			}
		}
	}
	else
	{	
		if(header)
		{
			li<-grep(end,the_names,ignore.case=T)
			if (length(li)>=1)
				end<-li[1];
		}
		if (is.na(end<-as.integer(end)))
			stop("The end column number is not an integer number nor a header name.")
		else if (end < 1)
			stop("The end column number is not a natural number.")
	}

	if (!is.na(end)) columns_accounted[end]=T

	width_names<-c('Width','Length')

	if (is.na(width))
	{
		if (header)
		{
			for (name in width_names)
			{
				li<-grep(name,the_names,ignore.case=T)
				if (length(li)<1) next
				width<-li[1];
				break;
			}
		}
	}
	else
	{	
		if(header)
		{
			li<-grep(width,the_names,ignore.case=T)
			if (length(li)>=1)
				end<-li[1];
		}
		if (is.na(width<-as.integer(width)))
			stop("The width column number is not an integer number not a header name.")
		else if (width < 1)
			stop("The width column number is not a natural number.")
	}

	if (!is.na(width)) columns_accounted[width]=T

	if (is.na(end) && is.na(width))
		stop("Both end and width parameters are not defined.")

	if (!is.na(end) && !is.na(width))
	{
		warning("Both end and width parameters are defined; width is ignored.")
		width=NA
	}

	arglist=(list(...))

	ctrl_names<-c('space','start','end','width')
	read_array<-c(space,start,end,width) # it is a list for futher checks
	names(read_array)<-ctrl_names


	for (a in names(arglist))
	{
		if (!is.na(arglist[[a]]))
		{
			col=arglist[[a]]
			if(header)
			{
				li<-grep(col,the_names,ignore.case=T)
				if (length(li)>=1)
					col<-li[1];
			}
			if (is.na(col<-as.integer(col)))
			{
				warning(paste("The column number for ",a," is not an integer number. It is ignored.",sep=""))
				next
			}
			else if (col < 1)
			{
				warning(paste("The column number for ",a," is not a natural number. It is ignored",sep=""))
				next
			}
			arglist[a]<-col
			read_array <- c(read_array,col)
			names(read_array)[length(read_array)] <- a
			columns_accounted[col]=TRUE
		}
		else
		{
				warning(paste("The column number for ",a," is NA. It is ignored.",sep=""))
				next
		}
	}

	if (header && autofeatures)
	{
		auto_added_features_pos<-(1:length(the_names))[!columns_accounted]
		for (i in auto_added_features_pos)
		{
			arglist=c(arglist,i)
			names(arglist)[length(arglist)]=the_names[i]
			read_array <- c(read_array,i)
			names(read_array)[length(read_array)] <- the_names[i]
		}

	}

	sorted_read_array=read_array[!is.na(read_array)] #remove na
	sorted_read_array=sorted_read_array[order(sorted_read_array)] #sort
	for (i in 1:(length(sorted_read_array)-1))
	{
		if (sorted_read_array[i]==sorted_read_array[i+1])
			warning(paste("The column numbers for ",names(sorted_read_array)[i]," and ",names(sorted_read_array)[i+1]," are equal. Are you sure?",sep=""))
	}

	#Now, we know that we have all necesssary parameters and 
	#all the parameters we have are natural and 
	#they are different (otherwise, the user is warned)


	ignored_list=c()

	for (a in names(read_array))
	{
		if (is.na(read_array[a]))
		{
			ignored_list=c(ignored_list,a)
			next
		}
		col=read_array[a]
		if (col>num_col)
		{
			if (length(ctrl_names[a==ctrl_names]>0) || ((is.na(end) && a=='width'))) # a is one of ctrl_names
				stop(paste("The column number for ",a," (",col,") is more than the number of columns.",sep=""))
			warning(paste("The column number for ",a," (",col,") is more than the number of columns. It is ignored.",sep=""))
			ignored_list=c(ignored_list,a)
			next
		}

	}

	mask=rep(TRUE,length(read_array))
	mask[ignored_list]=FALSE
	read_array=read_array[mask]

	#now, we removed all the oversized parameters

	#whether there are lines to skip (header)
	to_skip=TRUE
	skipped=0

	if (!is.na(end))
		second=end
	else
		second=width

	if(!is.integer(q[,start]))
		stop(paste("The column ",start," is not integer; maybe you forgot to skip some lines?.",sep=""))
	if(!is.integer(q[,second]))
		stop(paste("The column ",second," is not integer; maybe you forgot to skip some lines?.",sep=""))

	if (!is.na(end)) 
		ira=(IRanges(start=q[,start],end=q[,end]))
	else 
		ira=(IRanges(start=q[,start],width=q[,width]));

	#we created Iranges

	if (is.na(space) && length(list(...))==0) #it is an IRange
		return (ira);
		
	#if we are here, it is an RangedData
	if (is.na(space)) # the RangedData is for additional varaables; the name od the space is the name of the file
		rdata=RangedData(space=file,ira)
	else  # RangedData has space
		rdata=RangedData(space=q[,space],ira)
	for (a in names(arglist))
		rdata[[a]]=q[,arglist[[a]]]
	return(rdata)
}




