default: user_guide.pdf lpt.pdf cosmology_equations.pdf user_guide_public.pdf

%_public.pdf: %.tex
	pdflatex -jobname=$(basename $@) '\def\public{} \input{$<}'
	pdflatex -jobname=$(basename $@) '\def\public{} \input{$<}'

%.pdf: %.tex
	pdflatex $<
	pdflatex $<

clean :
	rm -f *.dvi *.aux *.bbl *.blg *.log

