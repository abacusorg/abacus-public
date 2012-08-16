default: cosmology_equations.pdf lpt.pdf

%.pdf: %.ps
	ps2pdf -dMaxSubsetPct=100 -dCompatibilityLevel=1.2 -dSubsetFonts=true -dEmbedAllFonts=true $<

%.ps: %.dvi
	dvips -Ppdf -t letter $< -o

%.dvi: %.tex
	latex $<
	latex $<
	/bin/rm $*.log

