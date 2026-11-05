SPHINXBUILD ?= sphinx-build
SPHINXOPTS ?=
DOXYGEN ?= doxygen

SOURCEDIR = source
BUILDDIR = _build
DOXYFILE = Doxyfile

.PHONY: all html doxygen clean help

all: html

html: doxygen
	$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" -W --keep-going $(SPHINXOPTS)

doxygen:
	mkdir -p "$(BUILDDIR)"
	$(DOXYGEN) "$(DOXYFILE)"

clean:
	rm -rf "$(BUILDDIR)" "$(SOURCEDIR)/api/generated"

help:
	@echo "make html     Build Doxygen XML and the Sphinx HTML site"
	@echo "make doxygen  Build only the Doxygen XML output"
	@echo "make clean    Remove generated documentation output"
