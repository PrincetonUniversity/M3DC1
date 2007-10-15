TEX = doc
BIB = doc
PDF = doc

ifdef BIB
	DEPS = ${TEX}.tex ${TEX}.bbl
else
	DEPS = ${TEX}.tex
endif

all : ${PDF}.pdf

${TEX}.bbl : ${BIB}.bib ${TEX}.aux
	bibtex ${TEX}
	latex ${TEX}

${TEX}.aux :
${TEX}.dvi : ${DEPS}
	latex $<

${TEX}.ps : ${TEX}.dvi
	dvips $< -o $@ -Ppdf

${PDF}.pdf : ${TEX}.ps
	ps2pdf $<

