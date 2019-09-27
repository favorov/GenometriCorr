# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2019 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#
# permuteIRanges is a function that permute randomly an IRanges


.permuteIRanges<-function(ir,chrom.length=NA)
{
	if (is.na(chrom.length)) chrom.length=max(end(ir))
	else if (chrom.length<max(end(ir)))
		stop("The IRanges given for permutation has a range in the that sticks out given chromosome lentgh.
  Not sure what to do.\n")		
	widths<-width(ir)
	result<-IRanges(start=sample(chrom.length, length(ir), replace=TRUE),width=widths) #permutting
	#from here until sort and return, we just care to cyclify the ranges that sticks out from chrom.length
	overlist<-end(result)>chrom.length 
	add_res<-IRanges()
	#print(result)
	for(index in which(overlist))
	{
		#cat(index,"\n")
		lag<-end(result)[index]-chrom.length
		end(result)[index]=chrom.length
		add_res<-c(add_res,IRanges(start=c(1),width=c(lag)))
	}
	result<-c(result,add_res)
	return(result[order(result)])
}


.rearrangeIRanges<-function(ir,chrom.length=NA)
{
	right_end=max(end(ir))
	if (is.na(chrom.length)) chrom.length=right_end
	else if (right_end>chrom.length)
		stop("The IRanges given for random rearramgement has a range in the that sticks out given chromosome lentgh.
  Not sure what to do.\n")		
	if(max(coverage(ir))>1)
		stop("You asked to rearrange an IRanges with overlaps.\n  Not sure what to do.\n")
	widths<-width(ir)
	addgap<-min(start(ir))+chrom.length-right_end #it is the sum of the left and right flanks
	gappes<-c(width(gaps(ir)),addgap) #we add it to the tail or the gaps
	widths<-widths[sample(length(widths))] #rearranging randomly
	gappes<-gappes[sample(length(gappes))] #rearranging randomly
	#print(widths)
	#print(gappes)
	starting_point=sample(tail(gappes,1),1)
	gappes=gappes[-length(gappes)]
	# the last gap is 'between' last and first ranges, so we sample starting point from it and stop using it
	result<-successiveIRanges(widths,gappes,starting_point)
	return (result)
}


