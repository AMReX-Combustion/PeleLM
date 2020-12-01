EPStoPDF = epstopdf

ALL: PeleLMUserGuide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = GettingStarted/PeleLMGettingStarted.tex \
	          Inputs/PeleLMInputs.tex \
	          Equations/PeleLMEquations.tex \
                  UnitsAndConstants/PeleLMUnits.tex \
                  AMR/AMR.tex \
                  Visualization/PeleLMVisualization.tex \
                  Software/PeleLMSoftware.tex \
                  Parallel/Parallel.tex \
                  Suggestions/Suggestions.tex \
                  Preface/PeleLMPreface.tex \
                  refs.bib

PeleLMUserGuide.dvi: PeleLMUserGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex PeleLMUserGuide.tex < /dev/null
	bibtex PeleLMUserGuide
	latex PeleLMUserGuide.tex < /dev/null
	latex PeleLMUserGuide.tex < /dev/null

PeleLMUserGuide.pdf: PeleLMUserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex PeleLMUserGuide.tex < /dev/null
	bibtex PeleLMUserGuide
	pdflatex PeleLMUserGuide.tex < /dev/null
	pdflatex PeleLMUserGuide.tex < /dev/null

pdf:	PeleLMUserGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) PeleLMUserGuide.dvi
	$(RM) PeleLMUserGuide.blg PeleLMUserGuide.log PeleLMUserGuide.toc
	$(RM) PeleLMUserGuide.aux PeleLMUserGuide.ps PeleLMUserGuide.bbl PeleLMUserGuide.pdf
	$(RM) PeleLMUserGuide.lot PeleLMUserGuide.lof PeleLMUserGuide.out
	$(RM) *~

.PHONY: clean
