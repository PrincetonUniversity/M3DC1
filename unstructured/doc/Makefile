TEX = doc
BIB = doc
PDF = doc

DEPS_COMMON = ${TEX}.tex SCORECdoc.tex
ifdef BIB
	DEPS = ${DEPS_COMMON} ${TEX}.bbl
else
	DEPS = ${DEPS_COMMON}
endif

all : ${PDF}.pdf

${TEX}.bbl : ${BIB}.bib ${TEX}.aux
	bibtex ${TEX}
	latex ${TEX}

${TEX}.aux : ${DEPS_COMMON}
	latex $<

${TEX}.dvi : ${DEPS}
	latex $<

${TEX}.ps : ${TEX}.dvi
	dvips $< -o $@ -Ppdf

${PDF}.pdf : ${TEX}.ps
	ps2pdf $<

clean :
	rm ${PDF}.pdf ${TEX}.ps ${TEX}.dvi ${TEX}.aux ${TEX}.bbl
