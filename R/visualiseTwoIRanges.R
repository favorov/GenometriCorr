# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2014 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors

# HJ -- pdf and close.device arguments

VisualiseTwoIRanges<-function(irA, irB, start=1, end=NA, nameA='RangesA', nameB='RangesB',
	chrom_length=NA, title=NA, pdf=NULL, close.device=NULL)
{
	pixelizeIRanges<-function(iranges,start,end,max.pixels=10000)
	{
		len=end-start+1
		if (len>max.pixels)
		{
			#less than max_pixels*bin are to cover len
			bin <- len %/% max.pixels	
			if ((len %% max.pixels) > 0) bin<-bin+1
			test.ranges<-IRanges(start=seq.int(from=1,by=bin,length.out=max.pixels),width=bin)
			overrle<-findOverlaps(test.ranges,iranges)
			masker<-tapply(width(iranges[subjectHits(overrle)]),queryHits(overrle),sum)
			mask<-rep(0,max.pixels)
			mask[as.integer(names(masker))]<-as.integer(masker)
		} else
		{
			mask<-as.vector(coverage(iranges,shift=1-start,width=len))
			bin<-1
		}
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
	# pdf
	if (!is.null(pdf)) {
		if (length(grep("\\.pdf$", pdf)) == 0)
			pdf <- paste(pdf, ".pdf", sep="")
		pdf(pdf)
		if (is.null(close.device)) close.device=TRUE
	} else
	{
		if (is.null(close.device)) close.device=FALSE
	}

	par(yaxt='n')
	plot(c(start,end), c(-1, 1), type = "n", xlab="", ylab="", bty='n')
  
  #find intersection for IRanges
  intersectionC<-intersect(irA, irB)

	#we are ready to think about what to plot.. \
	#let's get the length of the raster we are to prepare

	pixels<-as.integer(dev.size('px')[1]*(end-start+1)/xinch(dev.size('in')[1]))	
	#dev.size('px')[1] is how many pixels the thing will provide us with
	#end+start-1 is real (in bp) length of our plot
	#xinch is how many bp are there in argument inches
	#so, xinch(dev.size('in')[1])/dev.size('px')[1] 
	#is how many nucleotides are there in a pixel

	maskA<-pixelizeIRanges(irA,start,end,max.pixels=pixels)
	
	maskB<-pixelizeIRanges(irB,start,end,max.pixels=pixels)
  
	maskC<-pixelizeIRanges(intersectionC,start,end,max.pixels=pixels)
		
	#to debug	
	#maskA=sample(c(0,100,10000),pixels,replace=T,c(8/10,1/10,1/10))
	#maskB=sample(c(0,100,10000),pixels,replace=T,c(8/10,1/10,1/10))

	#
	#when we convert the data to pixels, we inversely regulate the coplementary colors.
	#So, (e.g. red channel) zero gives maximal level of the complement (see as.raster call below) 
	#and the color will be white. 
	#And, the maximal level of mask give zero to G and B channels, so we have pure red
	#For blue, we regulate G and R; 
	#For purple, only green is under regulation (Katya? is it correct?) 
 

	maxA<-max(maskA)
	if(maxA==0)
	{
		maxintenseA<-0
	} else
	{
		maskA <- 1-(maskA/maxA)
		#too sharp: maskA <- exp(-maskA) 
		maxintenseA<-max(maskA)
	}


	maxB<-max(maskB)
	if(maxB==0)
	{
		maxintenseB<-0
	} else
	{
		maskB <- 1-(maskB/maxB)
		#too sharp: maskB <- exp(-maskB) 
		maxintenseB<-max(maskB)
	}
  
	#maxC<-max(maskC)
	maxC<-min(c(maxA,maxB))
	if(maxC==0)
	{
	  maxintenseC<-0
	} else 
	{
	  maskC <- 1-(maskC/maxC) 
	  maxintenseC<-max(maskC)
	}
	
	img_len<-length(maskA)
	if (img_len != length(maskB))
	{
		stop("Image lengths are different, something went wrong.");
	}
	if (maxintenseA==0)  #not to divide by 0 in as.raster
		image_red<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_red<-
			as.raster(array(c(rep(maxintenseA,img_len),maskA,maskA),c(1,img_len,3)),max=maxintenseA)
  
  if (maxB==0)  #not to divide by 0 in as.raster
		image_blue<-
			as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
	else
		image_blue<-
			as.raster(array(c(maskB,maskB,rep(maxintenseB,img_len)),c(1,img_len,3)),max=maxintenseB)
  
  if (maxC==0) #not to divide by 0 in as.raster
    image_purple<-
	    as.raster(array(c(rep(1,img_len),rep(1,img_len),rep(1,img_len)),c(1,img_len,3)),max=1)
  else
    image_purple<-
	    as.raster(array(c(rep(1,img_len),maskC,rep(1,img_len)),c(1,img_len,3)),max=maxintenseC)
	
	#rasterImage(image_red, start, .05,end,.75)
	#rasterImage(image_blue, start, -.75,end,-.05)
	rasterImage(image_red, start, 0.20,end,0.75)
	rasterImage(image_blue, start, -0.75,end,-0.20)
	rasterImage(image_purple, start, 0.15,end,-0.15)
  
	text(c(start+len/2,start+len/2),c(.9,-.9),c(nameA,nameB))
	text(c(start-len/50),c(0,0),c('intersection'),srt=90)
	text(c(start+len+len/50),c(0.5),c(maxA),srt=90)
	text(c(start+len+len/50),c(0),c(maxC),srt=90)
	text(c(start+len+len/50),c(-0.5),c(maxB),srt=90)
	text(c(start+len/2),c(-1),c(paste0('Chromosome position: ',as.integer(xinch(dev.size('in')[1])/dev.size('px')[1]),' nucleotides per pixel')))
	if (!is.na(title))
	{
		title(main=title)
	}
	title(sub=paste0("On the right: number of covered nucleotides per pixel that yeilds full color intensity"))

	if (close.device)
		invisible(dev.off())
}
