PDF = a.pdf

only: ; cat a.tex | pdflatexfilter > a.pdf

rerun: $(PDF)

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	cat $< | pdflatexfilter > $@

clean : ; $(RM) $(PDF)
