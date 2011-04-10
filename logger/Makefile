PDFTEX = pdflatex

DOT = dot
#BIBTEX = bibtex

#TEX_BIBS = int.bib
#TEX_BBL := $(patsubst %.bib, %.bbl, $(TEX_BIBS))

# idea from http://ctan.unsw.edu.au/help/uk-tex-faq/Makefile
%.pdf: %.aux
	while ($(PDFTEX) $(patsubst %.aux, %.tex, $<) && \
		grep -q "Rerun to get cross" $(patsubst %.aux, %.log, $<)) do true ; \
	done

%.png: %.dia
	dia -t png -e $@ $<

%.aux: %.tex
	$(PDFTEX) $<

%.bbl: %.bib
	$(BIBTEX) $*

%.eps: %.dot
	dot -Teps $< > $@

.PHONY: all
all: logg5.pdf

logg5.pdf: logg5.tex Makefile

.PHONY: clean
clean:
	rm -rf $(patsubst %.tex, %.aux, $(TEX_SRCS))
	rm -rf $(patsubst %.tex, %.log, $(TEX_SRCS))
	rm -rf $(patsubst %.tex, %.out, $(TEX_SRCS))
	rm -rf $(patsubst %.tex, %.bbl, $(TEX_SRCS))
	rm -rf $(patsubst %.tex, %.blg, $(TEX_SRCS))

.PHONY: distclean
distclean: clean
	rm -rf logg5.pdf

