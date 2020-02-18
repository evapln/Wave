.PHONY: all clean help

FILENAME = rapport_wave

all:
	pdflatex $(FILENAME).tex
	pdflatex $(FILENAME).tex

clean:
	rm -f $(FILENAME).log $(FILENAME).aux $(FILENAME).toc *.gz

help:
	@echo "Usage:"
	@echo "  make [all]             Build the report"
	@echo "  make clean             Remove the generated files"
	@echo "  make help              Display the main targets of this Makefile with a short description"
