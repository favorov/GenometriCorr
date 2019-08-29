all: GenometriCorr.pdf clean 

GenometriCorr.pdf: \
	GenometriCorr.tex \
	png/*.png
	pdflatex GenometriCorr && pdflatex GenometriCorr && pdflatex GenometriCorr 

GenometriCorr.tex: GenometriCorr.Rnw
#ifdef CACHE
#	echo "library(weaver); Sweave('GenometriCorr.Rnw', driver=weaver())" | R --no-save --no-restore
#else
	echo "if ('unset'==Sys.getenv('CACHE',unset='unset')) {Sweave('GenometriCorr.Rnw')} else {library(weaver); Sweave('GenometriCorr.Rnw', driver=weaver())}" | R --no-save --no-restore
#endif

clean:
	rm -f *_picture.eps
	rm -f *_picture.pdf
	rm -f *_picture_vis.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.tex
	rm -f weaver_debug_log.txt
	rm -f *~


clean_all: clean
	rm -fr r_env_cache
