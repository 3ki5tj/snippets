prj=ocp

figs=
pdffigs=$(patsubst %,fig/%.pdf, $(figs))
epsfigs=$(patsubst %,fig/%.eps, $(figs))
figgps=$(patsubst %,fig/%.gp,  $(figs))

all: $(prj).pdf

$(prj).pdf: $(prj).tex liquid.bib Makefile $(pdffigs)
	#$(MAKE) -C fig
	pdflatex $<
	#bibtex $(prj)
	#pdflatex $<
	#pdflatex $<

$(pdffigs): %.pdf: %.gp
	$(MAKE) -C fig

clean:
	$(RM) *.aux *.bbl *.blg *.log *.dvi *.out \
	      *Notes.bib *~
	rstrip.py

$(prj)doc.zip: $(prj).tex $(prj).pdf liquid.bib \
	       $(epsfigs) $(pdffigs) $(figgps) \
	       README Makefile
	zip -r $@ $^

zip: $(prj)doc.zip

Bossman: zip
	rsync -avL $(prj)doc.zip \
	  cz1@129.109.88.204:/Bossman/cz1/vir/iedoc
