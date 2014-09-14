.PHONY: all clean try_clean_all clean_all 

cache="no"
#default

ifdef CACHE
cache="yes"
endif

# we have started from R
ifdef R_SESSION_TMPDIR 
do_clean_all="yes"
endif

ifdef R_USER 
do_clean_all="yes"
endif

ifdef R_HOME
R_CMD=$(R_HOME)/bin/R
else
R_CMD=R
endif

ifdef NOCLEAN
do_not_clean="yes"
endif

ifdef NO_CLEAN
do_not_clean=
endif

all: GenometriCorr.pdf try_clean_all

GenometriCorr.pdf: \
	GenometriCorr.tex \
	png/*.png
	pdflatex GenometriCorr && pdflatex GenometriCorr && pdflatex GenometriCorr 

GenometriCorr.tex: GenometriCorr.Rnw
ifeq ($(cache),"yes")
	echo "library(weaver); Sweave('GenometriCorr.Rnw', driver=weaver())" | $(R_CMD) --no-save --no-restore
else
	echo "Sweave('GenometriCorr.Rnw')" | $(R_CMD) --no-save --no-restore
endif

clean:
	rm -f *_picture.eps
	rm -f *_picture.pdf
	rm -f *_picture_vis.pdf
	rm -f *.aux
	rm -f *.log
	rm -f *.tex
	rm -f weaver_debug_log.txt
	rm -f *~
#	rm -f Rplots.pdf

try_clean_all:
ifndef do_not_clean
ifdef do_clean_all
	make clean_all
else
	make clean
endif
endif

clean_all: clean
	rm -fr data-examples
	rm -fr r_env_cache
