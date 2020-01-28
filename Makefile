.PHONY: all wave report clean help

all: wave report

wave:
	@cd implementation && $(MAKE)
	@cp -f implementation/wave ./

report:
	@cd rapport && $(MAKE)
	@cp -f rapport/rapport_wave.pdf ./

clean:
	@cd implementation && $(MAKE) clean
	@rm -f wave
	@cd rapport && $(MAKE) clean
	@rm -f rapport_wave.pdf

help:
	@echo "Usage :"
	@echo "  make [all]     Run the whole build of sudoku and report"
	@echo "  make wave      Build the executable file wave recursively in the directory implementation/"
	@echo "  make report    Build the report recursively in the directory rapport/"
	@echo "  make clean     Remove all files produced by the compilation"
	@echo "  make help      Display the main targets of this Makefile with a short description"
