EPStoPDF = epstopdf

ALL: LMCUserGuide.pdf 

eps_source = $(wildcard *.eps) $(wildcard */*.eps)

pdf_source = $(eps_source:.eps=.pdf)

section_sources = Introduction/LMCIntroduction.tex \
                  GettingStarted/LMCGettingStarted.tex \
	          Inputs/LMCInputs.tex \
	          Equations/LMCEquations.tex \
                  UnitsAndConstants/LMCUnits.tex \
                  AMR/AMR.tex \
                  Visualization/LMCVisualization.tex \
                  Software/LMCSoftware.tex \
                  Parallel/Parallel.tex \
                  Suggestions/Suggestions.tex

LMCUserGuide.dvi: LMCUserGuide.tex $(eps_source) $(section_sources) $(eps_source)
	latex LMCUserGuide.tex < /dev/null
	bibtex LMCUserGuide
	latex LMCUserGuide.tex < /dev/null
	latex LMCUserGuide.tex < /dev/null

LMCUserGuide.pdf: LMCUserGuide.tex $(pdf_source) $(section_sources) $(pdf_source)
	pdflatex LMCUserGuide.tex < /dev/null
	bibtex LMCUserGuide
	pdflatex LMCUserGuide.tex < /dev/null
	pdflatex LMCUserGuide.tex < /dev/null

pdf:	LMCUserGuide.pdf 

%.ps: %.dvi
	dvips -t letter -o $@ $<

%.pdf: %.eps
	$(EPStoPDF) $<

clean:
	$(RM) $(pdf_source) LMCUserGuide.dvi
	$(RM) LMCUserGuide.blg LMCUserGuide.log LMCUserGuide.toc
	$(RM) LMCUserGuide.aux LMCUserGuide.ps LMCUserGuide.bbl LMCUserGuide.pdf
	$(RM) LMCUserGuide.lot LMCUserGuide.lof LMCUsetGuide.out
	$(RM) *~

.PHONY: clean
