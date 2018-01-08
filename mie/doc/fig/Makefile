GP = $(wildcard *.gp)
PNGGP =
EPSGP = $(filter-out $(PNGGP), $(GP))
EPS = $(patsubst %.gp,%.eps, $(EPSGP))
PDF = $(patsubst %.gp,%.pdf, $(EPSGP))
PNG = $(patsubst %.gp,%.png, $(PNGGP))
XPDF =
epsdeps = Makefile $(wildcard ../data/*/*.dat)
pngdeps = Makefile

all: $(PNG) $(PDF) $(XPDF)

$(PNG) : %.png : %.gp $(pngdeps)
	gnuplot $<

$(PDF) : %.pdf : %.eps $(epsdeps)
	epstopdf $<

$(EPS) : %.eps : %.gp $(epsdeps)
	gnuplot $<

$(XPDF): figs.tex
	pdflatex -shell-escape $<
	mv figs-figure0.pdf massq.pdf
	mv figs-figure1.pdf vardep.pdf

clean:
	rm -f $(PDF) $(EPS) $(PNG) *~
	rstrip.py -Rv

