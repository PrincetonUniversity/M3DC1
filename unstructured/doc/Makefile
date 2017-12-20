TEX = doc
BIB = doc
PDF = doc

DEPS_COMMON = ${TEX}.tex SCORECdoc.tex idl_postprocessor.tex
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

${PDF}.pdf : ${TEX}.dvi
	dvipdf $<

clean :
	rm -f ${PDF}.pdf ${TEX}.dvi ${TEX}.aux ${TEX}.bbl
