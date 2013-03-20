default: user_guide.pdf lpt.pdf cosmology_equations.pdf

%.pdf: %.tex
	pdflatex $<
	pdflatex $<


clean :
	rm -f *.dvi *.aux *.bbl *.blg *.log

