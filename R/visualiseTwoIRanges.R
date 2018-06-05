# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2018 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
# VisualiseTwoIRanges is a function that visualise a pair of IRanges on a chromosome in two colors

# HJ -- pdf and close.device arguments


#'VisualiseTwoIRanges
#'
#'VisualiseTwoIRanges is a function that displays the intervals of two IRanges, one above the other, in different colors, along a chromosome or subset of a chromosome. The intent is to show large-scale relationships between the two IRanges.
#'
#'@param irA First \code{IRanges} object to be vsualised.
#'@param irB Second \code{IRanges} object to be vsualised.
#'@param start Start of the visualisation band.
#'@param end End of the visualisation band.
#'@param nameA Name of \code{irA}, default is 'RangesA'
#'@param nameB Name of \code{irB}, default is 'RangesB'
#'@param chrom_length The length of the chromosome spanned by \code{irA} and \code{irB}
#'@param title Title, printed at the top of the plot.
#'@param pdf Name of a file to which the image should be written. If \code{pdf=""} the filename is constructed from \code{nameA} and \code{nameB}. The suffix ".pdf" will be appended if not included. If \code{NULL}, no pdf is opened and an x11 window is raised up as by \code{plot}. 
#'@param close.device Whether to close the plot device after writing the image. A \code{FALSE} setting allows multiple images to be written to the same pdf file, but then eventually closing the device is up to the user. Default is \code{NULL} that means: close an x11, do not close a pdf.
#'@author Alexander Favorov \email{favorov@@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#'@references \href{http://genometricorr.sourceforge.net/}{GenometriCorr home}
#'@seealso The \code{\link{GenometriCorr}} documentation and vignette.
#'@examples
#'
#' library('rtracklayer')
#' library('GenometriCorr')
#' 
#' cpgis<-as(import(system.file("extdata", "UCSCcpgis_hg19.bed", package = "GenometriCorr")),'RangedData');
#' 
#' refseq<-as(import(system.file("extdata", "UCSCrefseqgenes_hg19.bed", package = "GenometriCorr")),'RangedData');
#' 
#' 
#' human.chrom.length<-c(249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,59373566,155270560)
#' 
#' 
#' names(human.chrom.length)<-c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrY','chrX')
#' 
#' VisualiseTwoIRanges(cpgis['chr1']$ranges, refseq['chr1']$ranges,
#' 	nameA='CpG Islands', nameB='RefSeq Genes',
#' 	chrom_length=human.chrom.length[['chr1']],
#' 	title="CpGIslands and RefGenes on chr1, hg19",
#' 	pdf='CpGi_vs_RefSeq_genes_chr1_hg19')
#'
#'@family GenometriCorr 2-range visualisation
#'@keywords hplot


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
			masker<-tapply(width(iranges[S4Vectors::subjectHits(overrle)]),S4Vectors::queryHits(overrle),sum)
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
		stop("The first parameter for TwoIRangesIndependence is to be IRanges.")
	
	if (!inherits(irB,"IRanges"))
		stop("The second parameter for TwoIRangesIndependence is to be IRanges.")

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
		if(pdf==".pdf")
			pdf<-paste0(nameA,"_",nameB,pdf)
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
	title(sub=paste0("On the right: covered nucleotides that yields full color intensity of a pixel"))

	if (close.device)
		invisible(dev.off())
}


#'VisualiseTwoGRanges
#'
#'VisualiseTwoGRanges is a function that displays the markups that are contatined two GRanges, in a multi-page pdf, page per chromosome, by \code{\link{VisualiseTwoIRanges}} for aech chromosome. The intent is to show large-scale relationships between the two markups. All the chromosome lenght data is provided by the GRanges.
#'
#'@param grA First \code{GRanges} object to be vsualised.
#'@param grB Second \code{GRanges} object to be vsualised.
#'@param nameA Name of \code{grA}, default is 'RangesA'
#'@param nameB Name of \code{grB}, default is 'RangesB'
#'@param title Title, printed at the top of the plot.
#'@param pdf Name of a file to which the image should be written. If \code{pdf=""} the filename is constructed from \code{nameA} and \code{nameB}. The suffix ".pdf" will be appended if not included. If \code{NULL}, no pdf is opened and an x11 window is raised up as by \code{plot}. 
#'@param close.device Whether to close the plot device after writing the image. A \code{FALSE} setting allows multiple images to be written to the same pdf file, but then eventually closing the device is up to the user. Default is \code{NULL} that means: close an x11, do not close a pdf.
#'@author Alexander Favorov \email{favorov@@sensi.org}, Loris Mularoni, Yulia Medvedeva, Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#'@references \href{http://genometricorr.sourceforge.net/}{GenometriCorr home}
#'@seealso The \code{\link{GenometriCorr}} documentation and vignette.
#'@family GenometriCorr 2-range visualisations
#'@keywords hplot

VisualiseTwoGRanges<-function(grA, grB, nameA='RangesA', nameB='RangesB', title=NA, pdf=NULL, close.device=NULL)
{
	if (!inherits(grA,"GRanges"))
		stop("The first parameter for TwoGRangesIndependence is to be GRanges.")
	
	if (!inherits(grB,"GRanges"))
		stop("The second parameter for TwoGRangesIndependence is to be GRanges.")

	sqinf<-seqinfo(grA)
	if(!all.equal(seqinfo(grB),sqinf))
		stop("The seqinfo() grA and grB is supposed to return exactly the same.")
	
	# pdf
	if (!is.null(pdf)) {
		if (length(grep("\\.pdf$", pdf)) == 0)
			pdf <- paste(pdf, ".pdf", sep="")
		if(pdf==".pdf")
			pdf<-paste0(nameA,"_",nameB,pdf)
		pdf(pdf)
		if (is.null(close.device)) close.device=TRUE
	} else
	{
		if (is.null(close.device)) close.device=FALSE
	}
	for (chn in 1:length(sqinf))
		VisualiseTwoIRanges(
			ranges(grA[start(seqnames(grA))[chn]:end(seqnames(grA))[chn]]),
			ranges(grB[start(seqnames(grB))[chn]:end(seqnames(grB))[chn]]),
			chrom_length = seqlengths(sqinf)[chn],
			nameA=nameA,
			nameB=nameB,
			pdf=NULL,
			title=ifelse(is.na(title),seqnames(sqinf)[chn],paste0(title,", ",seqnames(sqinf)[chn]))
		)
	if (close.device) 
		invisible(dev.off())
}
