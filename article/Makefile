#PDF  = $(shell ls *.tex | sed 's/\.tex/.pdf/')
PDF = a.pdf

pdf: $(PDF)

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	cat $< | pdflatexfilter > $@

clean : ; $(RM) $(PDF)
