NAME = spimap-manual
TEX = $(NAME).tex
PDF = $(NAME).pdf
DVI = $(NAME).dvi


IMAGES = 


all: $(PDF)

paper:
	pdflatex $(TEX)

$(PDF): $(TEX)
	pdflatex $(TEX)

bib:
	pdflatex $(TEX)
	bibtex $(NAME)
	pdflatex $(TEX)
	pdflatex $(TEX)

cpbib:
	cp ~/papers/papers.bib $(NAME).bib

clean:
	rm -rf $(PDF) 

