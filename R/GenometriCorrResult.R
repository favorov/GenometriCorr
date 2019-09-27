# GenometriCorrelation project evaluating two interval markups genomewide independence. 
# (c) 2010-2019 Alexander Favorov, Loris Mularoni, Yulia Medvedeva, 
#               Harris A. Jaffee, Ekaterina V. Zhuravleva, Leslie M. Cope, 
#               Andrey A. Mironov, Vsevolod J. Makeev, Sarah J. Wheelan.
#
# result - result of GenometriCorrelation function - works with ini file format  

#if (!require('methods')) stop('GenometriCorrResult requires methods package!\n')
#if (!require('graphics')) stop('GenometriCorrResult requires graphics package!\n')

#'@importFrom grDevices as.raster rainbow
#'@importFrom grDevices colorRamp dev.off dev.size pdf rgb 
#'@importFrom graphics axis hist layout lines mtext par plot plot.new rasterImage text xinch 

#'@export
setClass('GenometriCorrResult',contains='namedList',representation(config="GenometriCorrConfig"))

#'@export
setMethod('show','GenometriCorrResult',function(object)
	{
		if (length(object) == 0) {
			cat("It is an empty GenometriCorrResult.")
			return()
		}
		namelist<-names(object[[1]])
		#possible operations with namelist here
		do_not_show<-
			c(
				'relative.distances.data',
				'absolute.min.distance.data',
				'absolute.inter.reference.distance.data',
				'relative.distances.ecdf.deviation.area.null.list',
				'scaled.absolute.min.distance.sum.null.list',
				'jaccard.measure.null.list',
				'jaccard.intersection.null.list',
				'projection.test')
				#scaled.absolute.min.distance.sum
		#we remove the do_not_show list from namelist
		namelist<-setdiff(namelist,do_not_show)
		print(sapply(object,function(x){
		x[namelist]
		}))
	})


#setGeneric('plot', useAsDefault=function(x,y,...) graphics::plot(x,y,...))
setGeneric('graphical.report',function(x,pdffile='',show.all=TRUE,show.chromosomes=c(), trustname=TRUE, make.new=TRUE)
                                standardGeneric('graphical.report'))


#'@export
setMethod('graphical.report',
	signature(x='GenometriCorrResult'),
	function(x, pdffile='',show.all=FALSE,show.chromosomes=c(),trustname=TRUE,make.new=TRUE)
	# if show.all == TRUE, show all of them
	# if show.all == FALSE, show only awhole+show.chromosomes; if it leads to showing nothing, show.all turns back to TRUE
	{
		if (!is.null(x@config$options$awhole.space.name))
			awhole.space.name=x@config$options$awhole.space.name
		else
			awhole.space.name='awhole'


		if (is.null(x[[awhole.space.name]]))
			awhole.space.name.list=c()
		else
			awhole.space.name.list=c(awhole.space.name)	

		if(!show.all)
		{
			to.show<-union(show.chromosomes,awhole.space.name.list)
			if (length(to.show)==0) show.all=TRUE
		}

		if (show.all)
			to.show<-names(x)	

		chromosomes<-setdiff(names(x),awhole.space.name.list)

		show.awhole<-(length(awhole.space.name.list)>0)
		#we believe that both reference and exist, so we do not check it. 
		if (pdffile=='') 
		{
			if (!is.null(x@config$data$query) && length(grep("[a-zA-Z0-9]",x@config$data$query))>0)
			{
				if (length(strsplit(basename(x@config$data$query),'\\.')[[1]])>1)
					qname<-paste(head(strsplit(basename(x@config$data$query),'\\.')[[1]],-1),sep='.')
				else
					qname<-basename(x@config$data$query)
			}
			else qname<-'unknown'

			if (!is.null(x@config$data$reference) && length(grep("[a-zA-Z0-9]",x@config$data$reference))>0)
			{
				if (length(strsplit(basename(x@config$data$reference),'\\.')[[1]])>1)
					rname<-paste(head(strsplit(basename(x@config$data$reference),'\\.')[[1]],-1),sep='.')
				else
					rname<-basename(x@config$data$reference)
			}
			else rname<-'unknown'
			pdffile<-paste(qname,'_to_',rname,'_plot.pdf',sep='')
		}
		else if (trustname==FALSE)
		{
			if (!(length(strsplit(pdffile,'\\.')[[1]]) > 1 && tail(strsplit(pdffile,'\\.')[[1]],1) == 'pdf'))
				pdffile<-paste(pdffile,'.pdf',sep='')
		}
		


		width<-if(! is.null(x@config$options$keep.distributions) && x@config$options$keep.distributions) 12 else 4
		if (make.new == TRUE)
		{
			pdf(file=pdffile, width=width, height=4, paper="special")
		}
		for (name in to.show)
		{
			if (! is.null(x@config$options$keep.distributions) && x@config$options$keep.distributions)
				layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), heights=c(1,1,1))
			plot.new()
			data <- x[[name]]#[[1]]
			if (name == "awhole")
				usename <- "All Chromosomes"
			else
				usename <- name
			formstring_s<-"%s"
			formstring_g<-"%g"
			mtext(usename, line=-1, cex=1.5)
			mtext(paste("Query population :", data$query.population, sep=" "), line=-4, cex=0.7)
			mtext(paste("Reference population :", data$reference.population, sep=" "), line=-5, cex=0.7)
			mtext(paste("Relative Ks p-value :", sprintf(formstring_g,data$relative.distances.ks.p.value), sep=" "), line=-6, cex=0.7)
			mtext(paste("Relative ecdf deviation area :", sprintf(formstring_g,data$relative.distances.ecdf.deviation.area), sep=" "), line=-7, cex=0.7)
			mtext(paste("Relative ecdf area correlation :", sprintf(formstring_g,data$relative.distances.ecdf.area.correlation), sep=" "), line=-8, cex=0.7)
			if('relative.distances.ecdf.deviation.area.p.value' %in% names(data))
				mtext(paste("Relative ecdf deviation area p-value :", sprintf(formstring_s,data$relative.distances.ecdf.deviation.area.p.value), sep=" "), line=-9, cex=0.7)
			if('scaled.absolute.min.distance.sum.p.value' %in% names(data))
				mtext(paste("Scaled Absolute min. distance p-value :", sprintf(formstring_s,data$scaled.absolute.min.distance.sum.p.value), sep=" "), line=-10, cex=0.7)
				mtext(paste("Scaled Absolute min. lower tail :", data$scaled.absolute.min.distance.sum.lower.tail, sep=" "), line=-11, cex=0.7)
			if('jaccard.measure.p.value'%in% names(data))
				mtext(paste("Jaccard Measure p-value :", sprintf(formstring_s,data$jaccard.measure.p.value), sep=" "), line=-12, cex=0.7)
			if('jaccard.measure.lower.tail'%in% names(data))
				mtext(paste("Jaccard Measure lower tail :", data$jaccard.measure.lower.tail, sep=" "), line=-13, cex=0.7)
			if('projection.test.p.value'%in% names(data))
				mtext(paste("Projection test p-value :", sprintf(formstring_g,data$projection.test.p.value), sep=" "), line=-14, cex=0.7)
			if('projection.test.lower.tail'%in% names(data))
				mtext(paste("Projection test lower tail :", data$projection.test.lower.tail, sep=" "), line=-15, cex=0.7)
			if('projection.test.obs.to.exp'%in% names(data))
				mtext(paste("Projection test observed to expected ratio :", sprintf(formstring_g,data$projection.test.obs.to.exp), sep=" "), line=-16, cex=0.7)

			if (! is.null(x@config$options$keep.distributions) && x@config$options$keep.distributions)
			{
				twotimes<-function(x){2*x}
				plot(ecdf(data$absolute.min.distance.data),col="black",lty='solid',lwd=2, main="Absolute distances", xlab="Distance (bp)", ylab="Cumulative fraction")
				#plot(ecdf(data$absolute.inter.reference.distance.data), col="grey",lty='solid',lwd=2, main="Absolute distances", xlab="Distance (bp)", ylab="Cumulative fraction")
				#plot abs.data? no sense :)
				#create the null-hypothesis cdf function for absolute.min.distance.data
				cdf.x<-sort(data$absolute.inter.reference.distance.data)/2
				cdf.y<-cdf.x #to init
				L.half<-sum(cdf.x) #just to remember it is half-sum 
				for(i in 1:length(cdf.x)) #it is sorted, a=x_i=r/2; after the loop, we will 2/L
					cdf.y[i]<-sum(cdf.x[1:i])+ # x is r/2 and sum is L/2
						(length(cdf.x)-i)*cdf.x[i]
				cdf.y<-cdf.y/L.half
				lines(cdf.x,cdf.y,col="blue") #plot expetation with blue
				#ok. start plot 2 with expetation
				plot(twotimes, xlim=c(0,0.5), col="blue", main="Relative distances", xlab="Fractional distance", ylab="Cumulative fraction")
				lines(ecdf(data$relative.distances.data))
			}
		}
		if (make.new == TRUE)
		{
			dev.off()
		}
	})



setGeneric('visualize',function(x,pdffile='',show.all=TRUE,show.chromosomes=c(), trustname=TRUE, make.new=TRUE, style="blue-white-red")
                                standardGeneric('visualize'))

#'@export
setMethod('visualize',
	signature(x='GenometriCorrResult'),
	function(x, pdffile='',show.all=FALSE,show.chromosomes=c(), trustname=TRUE, make.new=TRUE, style="blue-white-red")
	{
		# gplots is imported now
		#if (!require('gplots')) stop('GenometriCorr visualize requires gplots package\n')
		if (is.null(x@config$options$keep.distributions)) return()

		if (!is.null(x@config$options$awhole.space.name))
			awhole.space.name=x@config$options$awhole.space.name
		else
			awhole.space.name='awhole'

		if (is.null(x[[awhole.space.name]]))
			awhole.space.name.list=c()
		else
			awhole.space.name.list=c(awhole.space.name)	

		if(!show.all)
		{
			to.show<-union(show.chromosomes,awhole.space.name)
			if (length(to.show)==0) show.all=TRUE
		}

		if (show.all)
			to.show<-names(x)

		chromosomes<-setdiff(names(x),awhole.space.name.list)

		show.awhole<-(length(awhole.space.name.list)>0)
		alldist <- c() # maybe, we wiil never need it
		alldatdist <- c()

		if(show.all == TRUE)
		{
			for (chr in chromosomes)
			{
				data <- x[[chr]]#[[1]]
				alldist <- c(alldist, diff(data$reference.middles))
				alldatdist <- c(alldatdist, data$absolute.min.distance.data)
			}
			alldist <- as.numeric(alldist)
			alldatdist <- as.numeric(alldatdist)
		}
		if (pdffile=='') 
		{
			if (!is.null(x@config$data$query))
			{
				if (length(strsplit(basename(x@config$data$query),'\\.')[[1]])>1)
					qname<-paste(head(strsplit(basename(x@config$data$query),'\\.')[[1]],-1),sep='.')
				else
					qname<-basename(x@config$data$query)
			}
			else qname<-'unknown'

			if (!is.null(x@config$data$reference))
			{
				if (length(strsplit(basename(x@config$data$reference),'\\.')[[1]])>1)
					rname<-paste(head(strsplit(basename(x@config$data$reference),'\\.')[[1]],-1),sep='.')
				else
					rname<-basename(x@config$data$reference)
			}
			else rname<-'unknown'
			pdffile<-paste(qname,'_to_',rname,'_visualize.pdf',sep='')
		}
		else if (trustname==FALSE)
		{
			if (!(length(strsplit(pdffile,'\\.')[[1]]) > 1 && tail(strsplit(pdffile,'\\.')[[1]],1) == 'pdf'))
				pdffile<-paste(pdffile,'.pdf',sep='')
		}

		if (make.new == TRUE)
		{
			pdf(file=pdffile, width=10, height=19, paper="special")
			mymat <- matrix(ncol=2, nrow=8)
			mymat[1,1] <- 2
			mymat[1,2] <- 2
			mymat[2,1] <- 3
			mymat[2,2] <- 3
			mymat[3,1] <- 1
			mymat[3,2] <- 1
			mymat[4,1] <- 4
			mymat[4,2] <- 5
			mymat[5,1] <- 6
			mymat[5,2] <- 6
			mymat[6,1] <- 7
			mymat[6,2] <- 8 
			mymat[7,1] <- 9 
			mymat[7,2] <- 9 
			mymat[8,1] <- 10
			mymat[8,2] <- 11
				
			layout(mymat, heights=c(0.05,0.05,0.15,0.15,0.15,0.15,0.15,0.15))
		}
		nbreaks <- 50
		colorinterpolation <- 10
		layoutresults <- 3
		i <- 0
		for (name in to.show)
		{
			maxdist <- 0
			i <- i+1
			data <- x[[name]]#[[1]]
			plot(1,1, type="n", axes=F, xlab="", ylab="")
			par(mar=c(1,1,1,1))
			if (name == "awhole")
			{
				name <- "All chromosomes"
			}
			mtext(paste("Results: ", name, sep=""), line=-2, cex=1.0, side=3)
			mtext(paste("Overlap summary (Jaccard and projection tests)"), line=-4, cex=0.8, side=3)
			mtext(paste("Jaccard p-value: ", data$jaccard.measure.p.value, sep=""), line=-6, cex=0.8, side=3)
			if (!is.integer(data$jaccard.measure.p.value))
			{
				if (data$jaccard.measure.lower.tail == FALSE)
				{
					mtext(paste("Query and reference intervals overlap significantly more than expected by chance, by Jaccard"), line=-8, cex=0.8, side=3)
				}
				else
				{
					mtext(paste("Query and reference intervals overlap significantly less than expected by chance, by Jaccard"), line=-8, cex=0.8, side=3)
				}
			}
			if (!is.integer(data$projection.test.p.value))
			{
				if (data$projection.test.lower.tail == FALSE)
				{
					mtext(paste("Query midpoints and reference intervals overlap significantly more than expected by chance, by projection"), line=-10, cex=0.8, side=3)
				}
				else
				{
					mtext(paste("Query midpoints and reference intervals overlap significantly less than expected by chance, by projection"), line=-10, cex=0.8, side=3)
				}
			}

			if (!is.integer(data$relative.distances.ecdf.area.correlation))
			{
				goodness_of_corr <- 0
			} else
			{
				goodness_of_corr <- data$relative.distances.ecdf.area.correlation
			}
			data_rel_exp <- seq(0, 50, by=50/nbreaks)/100
			data$relative.distances.data <- data$relative.distances.data[data$relative.distances.data>=0]
			hist_to_plot <- hist(data$relative.distances.data, breaks=data_rel_exp, plot=F)
			hist_exp <- hist(data_rel_exp, breaks=hist_to_plot$breaks, plot=F)
			plottitle <- paste("Relative distance from query to reference, pvalue= ", data$relative.distances.ecdf.deviation.area.p.value, sep="")
			.plot_colored_hist(goodness_of_corr, hist_to_plot$density, log2(hist_to_plot$density/hist_exp$density), nbreaks, colorinterpolation, plottitle, 0.5, printlegend = switch((i%%layoutresults)+1, 0, 1, 0, 0),style=style)
			if (!is.integer(data$scaled.absolute.min.distance.sum.p.value))
			{
				goodness_of_corr = 0
			} else
			{
				goodness_of_corr <- as.numeric(data$scaled.absolute.min.distance.sum.p.value)
			}
#for absolute distances the data are really different than for relative distances
#need to consider both the real interval lengths and the actual distances
			if (maxdist==0)
			{
				maxdist <- min(quantile(data$absolute.min.distance.data)[[4]], quantile(data$absolute.inter.reference.distance.data)[[4]])
				maxdist <- min(maxdist, min(quantile(data$absolute.min.distance.data)[[3]], quantile(data$absolute.inter.reference.distance.data)[[3]]))
			}
                        if (maxdist<nbreaks)
                        {
                                maxdist <- max(nbreaks,max(quantile(data$absolute.min.distance.data)[[4]], quantile(data$absolute.inter.reference.distance.data)[[4]]))
                        }
			dist_to_consider <- data$absolute.inter.reference.distance.data[data$absolute.inter.reference.distance.data<=maxdist]
			exp_dist_to_consider <- data$absolute.min.distance.data[data$absolute.min.distance.data<=maxdist]
			exp_dist_to_consider <- exp_dist_to_consider[exp_dist_to_consider >= 0]
			abs_breaks <- seq.int(from=0, to=maxdist+maxdist/nbreaks-1, by=maxdist/nbreaks)
			if (length(dist_to_consider) == 0 || length(exp_dist_to_consider) == 0)
			{
				plot(1, 1, type="n", axes=F, xlab="", ylab="")
				par(mar=c(1,1,1,1))
				mtext("Insufficient data", line=-3)
				next
			}
			hist_obs <- hist(exp_dist_to_consider, breaks=abs_breaks, plot=F)
			hist_exp <- hist(dist_to_consider, breaks=hist_obs$breaks, plot=F)
			abs_breaks <- hist_obs$breaks
			hist_exp$counts[hist_exp$counts<=0] <- min(hist_exp$counts[hist_exp$counts>0])
			hist_exp$counts <- hist_exp$counts*(sum(hist_obs$counts)/sum(hist_exp$counts))  #scale expected to same # observations
			plottitle <- paste("Absolute distance from query to reference, pvalue= ", data$scaled.absolute.min.distance.sum.p.value, sep="")
			.plot_colored_hist(goodness_of_corr, hist_obs$density, log2(hist_obs$counts/hist_exp$counts), length(abs_breaks), colorinterpolation, plottitle, maxdist, style=style)
		}
		if (make.new == TRUE)
		{
			dev.off()
		}
	})


.plot_colored_hist <- function(goodness_of_corr, densities, obs_exp_rats, nbreaks, colorinterpolation=1, plottitle, extremexlim, printlegend=0, style="blue-white-red")
{
	obs_exp_rats[!(is.finite(obs_exp_rats))] <- 0
	nbreaks=length(obs_exp_rats)
##first scale the correlation by the range of difference in strength of signal
        span_graph <- max(median(obs_exp_rats)-min(obs_exp_rats), max(obs_exp_rats)-median(obs_exp_rats))
        scale_factor_line <- 4/span_graph
        scale_factor <- 1-goodness_of_corr
	if (printlegend > 0)
	{
		par(mar=c(1,1,5,1))
		posscols_up <- .rainbooo(1000, start=0.05, end=0.18,style=style)
		posscols_down <- .rainbooo(1000, start=0.18, end=0.6,style=style)
		allposscols <- c(posscols_up, posscols_down)
		allposscols <- allposscols[length(allposscols):1]
		hist(seq(2:1999), breaks=(seq(1:2002)-0.5), col=allposscols, border=allposscols, axes=F, labels=F, main="Color key\n  <- blue is negative correlation, -> red is positive correlation", xlab="", ylab="", cex.main=1.2)
		par(mar=c(1,1,1,1))
		plot(1,1, type="n", axes=F, xlab="", ylab="")
		mtext(paste("Overlay line on graph is data density, over ", nbreaks, "  bins\nThis range of densities is real though does not on its own convey significance\nThe p-value signals whether the trends are statistically significant.", sep=""), line=-3, cex=1.0, side=3)
	}
        ranges_to_plot <- obs_exp_rats*scale_factor_line
        end <- 0.18+scale_factor*0.42
	#if (end > 1)
	#{
	#	end = 1
	#}
	if (end > .6)
	{
		end = .6
	} # consistent with posscols
	colorscheme_down <- .rainbooo(nbreaks, start=0.18, end=end,style=style) ##scaling by correlation strength
	ranges_to_plot_up <- ranges_to_plot[ranges_to_plot>=0]
	start=0.18-scale_factor*0.13
	#if (start < 0)
	#{
	#	start = 0.05
	#}
	if (start < 0.05)
	{
		start = 0.05
	}
	colorscheme_up <- .rainbooo(nbreaks, start=start, end=0.18,style=style) ##scaling by correlation strength
	##reverse the positive correlation colorscheme
	colorscheme_up <- colorscheme_up[nbreaks:1]
	ranges_to_plot_down <- ranges_to_plot[ranges_to_plot<0]
	thesecols_up <- 0
	thesecols_down <- 0
        if (length(ranges_to_plot_up) > 0)
        {
                rr <- max(ranges_to_plot_up)-min(ranges_to_plot_up)
                thesecols_up <- ((ranges_to_plot_up-min(ranges_to_plot_up))/rr)*(nbreaks-1)
                thesecols_up <- as.integer(thesecols_up+1)
        }
        if (length(ranges_to_plot_down) > 0)
        {
                rr <- max(ranges_to_plot_down)-min(ranges_to_plot_down)
                thesecols_down <- ((abs(ranges_to_plot_down)-min(abs(ranges_to_plot_down)))/rr)*(nbreaks-1)
                thesecols_down <- as.integer(thesecols_down+1)
        }
	mybreaks <- seq(1:(nbreaks+1))
	mybreaks <- mybreaks-0.5
	allcols <- vector(length=length(obs_exp_rats))
	allcols[ranges_to_plot<0] <- colorscheme_down[thesecols_down]
	allcols[ranges_to_plot>=0] <- colorscheme_up[thesecols_up]
	if (colorinterpolation > 0)
	{
		newcols <- vector(length=colorinterpolation*(length(allcols)-1))
		for (v in 1:(length(allcols)-1))
		{
			newcols[(v-1)*colorinterpolation + 1] = allcols[v]
			ramp <- gplots::colorpanel(colorinterpolation, allcols[v], allcols[v+1])
			for (j in 1:colorinterpolation)
			{
				newcols[(v-1)*colorinterpolation + 1 + j] = ramp[j]
			}
		}
		allcols <- newcols
		newbreaks <- vector(length=(length(mybreaks)-1)*colorinterpolation)
		newbreaks <- seq(1:(length(newbreaks)))
		mybreaks <- newbreaks-0.5
	}
	par(mar=c(3,3,1,1))
	hist(seq(1:(length(allcols)-1)), breaks=mybreaks[mybreaks<=(length(allcols))], col=allcols, border=allcols, axes=FALSE, labels=FALSE, main=plottitle, xlab="", ylab="", cex.main=1.2)
	par(new=TRUE)
	par(mar=c(3,3,1,1))
	uprange <- max(densities)/sum(densities)
	plot(seq(1:length(densities)), densities/sum(densities), ylim = c(0-uprange/10,uprange+uprange/10), type="l", xaxt="n", yaxt="n", ylab="", xlab="")
	xmax = length(obs_exp_rats)-1
	xlabs <- c(0, extremexlim/5, extremexlim*2/5, extremexlim*3/5, extremexlim*4/5, extremexlim)
	xlabs <- round(xlabs, 3)
	atlabs <- c(0, xmax/5, xmax*2/5, xmax*3/5, xmax*4/5, xmax)
	atlabs <- round(atlabs, 3)
	axis(1, labels=xlabs, at=atlabs)
	topdensity <- round(max(densities)/sum(densities), 3)
        topdensity_to_print <- topdensity*100
	tdlab <- paste(topdensity_to_print, "%", sep="")
	axis(2, labels=c(0,tdlab), at=c(0,topdensity))
}

.rainbooo<-function(n, start, end,  style="blue-white-red")
{
		
		#cat('.rainbbbooo    ',n,style,start,end,"\n")
		if(style=="rainbow")
		{
			return(rainbow(n=n,start=start,end=end))
		}
		if (style!="blue-white-red")
			warning("Unknown color style; use dafault blue-white-red")


		# remap any...0.18 to red (0.05->0) -> white(0.18->1)
		# remap 0.18..any to white(0.18->0) -> blue(0.6->1) 
		#3-color ramp code
		#ramp <- colorRamp(c('blue',"white",'red'))
		#frCol = rgb( ramp(seq(0, 1, length = 20)), max = 255)
		{
			if (end==0.18)
			{
				if (start<0.05)
				{
					warning("Color mapping from <0.05 to 0.18 in .rainbooo.\n")
					return(rainbow(n=n,start=start,end=end))
				}
				ramp <- colorRamp(c("red","white"))
				frCol = rgb( ramp(seq((start-0.05)/(0.18-0.05), 1, length = n)), maxColorValue = 255)
			}
			else if (start==0.18)
			{
				if (end>0.6)
				{
					warning("Color mapping from 0.18 to >0.6 in .rainbooo.\n")
					return(rainbow(n=n,start=start,end=end))
				}
				ramp <- colorRamp(c("white","blue"))
				frCol = rgb( ramp(seq(0, 1-((0.6-end)/(0.6-0.08)), length = n)), maxColorValue = 255)
			}
			else
			{
				warning("Color mapping not from 0.18 and not to 0.18 in .rainbooo.\n")
				return(rainbow(n=n,start=start,end=end))
			}
		}
}


