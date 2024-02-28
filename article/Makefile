PDF = a.pdf

only: ; cat a.tex | pdflatexfilter > a.pdf

rerun: a.tex
	grep ^%SCRIPT $< | cut -c9- | sh

%.pdf : %.tex
	grep ^%SCRIPT $< | cut -c9- | sh
	cat $< | bash pdflatexfilter > $@

clean : ; $(RM) $(PDF)
