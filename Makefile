PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

all: rd check clean

alldocs: rd readme


rd:
	Rscript -e 'roxygen2::roxygenise(".")'

readme:
	# Rscript -e 'rmarkdown::render("README.Rmd")'
	quarto render README.qmd

build: 
	#cd ..;\
	#R CMD build $(PKGSRC)
	Rscript -e 'devtools::build()'

build2: rd
	cd ..;\
	R CMD build --no-build-vignettes $(PKGSRC)

install: 
	#cd ..;\
	#R CMD INSTALL $(PKGNAME)_$(PKGVERS).tar.gz
	Rscript -e 'devtools::install()'

check: 
	# cd ..;\
	# Rscript -e 'rcmdcheck::rcmdcheck("$(PKGNAME)_$(PKGVERS).tar.gz", args="--as-cran")'
	Rscript -e 'devtools::check()'

check2: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz

bioccheck:
	cd ..;\
	Rscript -e 'BiocCheck::BiocCheck("$(PKGNAME)_$(PKGVERS).tar.gz")'

clean:
	cd ..;\
	$(RM) -r $(PKGNAME).Rcheck/

