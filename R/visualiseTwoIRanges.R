# GenometriCorrelation project evaluating two markups genomwide independence. 
# (c) 2010-2011 Alexander Favorov, Leslie Cope, Yulia Medvedeva, 
#              Loris Mularoni, Vsevolod Makeev, Sarah Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors
# $Id: visualiseTwoIRanges.R 1717 2012-04-18 19:29:15Z favorov $


VisualiseTwoIRanges<-function(irA,irB,start=1,end=NA, nameA='RangesA',nameB='RangesB',chrom_length=NA,title=NA)
{
	pixelizeIRange<-function(irange,start,end,max_pixels=10000)
	{
		len=end-start+1
		if (len>max_pixels)
		{	
			step=ceiling(len/max_pixels)
			seqpos=seq(start,end,step)
			#print(step)
			#print (seqpos)
			#slow way 1
			#mask<-rep(0,length(seqpos))
			#ipos=1
			#for (pos in seqpos)
			#{
			#	su<-sum(coverage(irange,shift=1-pos,width=step))
			#	mask[ipos]=su
			#	ipos=ipos+1
			#}
			#slow way2
			#mask<-sapply(
			#			seqpos,
			#			function(pos)
			#			{
			#				su=sum(coverage(irange,shift=1-pos,width=step));
			#				#cat(c(pos," ",step," ",su,"\n"));
			#				return (su)}
			#			)
			cover<-coverage(irange,shift=1-start,width=end)
			rW<-width(cover)
			rX<-start(cover)+start
			rE<-end(cover)+start
			N<-nrun(cover)
			rV<-runValue(cover)
			current<-1 #current run of RLE  object 'cover' 
			mask<-rep(0,length(seqpos))
			ipos=1
			for (pos in seqpos)
			{
				pixel_weight<-0
				while(rX[current]<=pos+step-1 && current<=N)
				{
					if (rE[current]>=pos) # they overlap
					{
						#just overlap of two ranges
						overlap.length <- min(rE[current],pos+step-1) - max(rX[current],pos)
						pixel_weight <- pixel_weight+overlap.length*rV[current]
					}
					current <- current+1
				}
				mask[ipos]=pixel_weight
				ipos=ipos+1
			}
		}
		else
		{
			mask<-as.vector(coverage(irange,shift=1-start,width=len))
			step<-1
		}
		#print (mask)
		return(mask)
	}

	#require(grDevices)

	if (!inherits(irA,"IRanges"))
		stop("The first parameter for TwoIRangesIndependence is to be IRange.")
	
	if (!inherits(irB,"IRanges"))
		stop("The second parameter for TwoIRangesIndependence is to be IRange.")

	if (is.na(chrom_length))
	{
		if (is.na(end)) 
		{
			chrom_length<-chromosomes.length.eval(irA,irB)
			end<-chrom_length
			warning("Chromosome length and the right margin of picture are evaluated rather than pre-given")
		}
		else
			chrom_length<-chromosomes.length.eval(irA,irB)
	}

	if (is.na(end))
		end<-chrom_length

	if (end<=start)
		stop("End is less or equal than start.")

	len<-end-start+1
	maskA<-pixelizeIRange(irA,start,end)
	
	#print(maskA)
	maskB<-pixelizeIRange(irB,start,end)
	
	#print(maskB)

	#maskA=sample(c(0,100,10000),10000,replace=T,c(8/10,1/10,1/10))
	#maskB=sample(c(0,100,10000),10000,replace=T,c(8/10,1/10,1/10))

	#maxA <- max(maskA)
	#maxB <- max(maskB)
	#maskA <- maxA-maskA
	#maskB <- maxB-maskB
	maskA <- exp(-maskA) #to be better on view, now 0 is white and everuthing > 0 is a bit or more red 
	maskB <- exp(-maskB) #to be better on view, now 0 is white and everuthing > 0 is a bit or more blue
 
	maxA<-max(maskA)
	maxB<-max(maskB)
	img_len<-length(maskA)
	if (img_len != length(maskB))
	{
		stop("Image lengthes are different, something went wrong.");
	}
	if (maxA==0)  #not to divide by 0 in as.raster
		image_red<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_red<-
			as.raster(array(c(rep(maxA,img_len),maskA,maskA),c(1,img_len,3)),max=maxA)
	if (maxB==0)  #not to divide by 0 in as.raster
		image_blue<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_blue<-
			as.raster(array(c(maskB,maskB,rep(maxB,img_len)),c(1,img_len,3)),max=maxB)
	par(yaxt='n')
	plot(c(start,end), c(-1, 1), type = "n", xlab="", ylab="")
	rasterImage(image_red, start, .05,end,.75)
	rasterImage(image_blue, start, -.75,end,-.05)
	text(c(start+len/2,start+len/2),c(.9,-.9),c(nameA,nameB))
	if (!is.na(title))
	{
		title(main=title)
	}
}
