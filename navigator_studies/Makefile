# Build the looper design note and its figures.
#
#   make        -> figures/loopers.pdf (the note, with figures)
#   make figs   -> the standalone TikZ figures only
#   make png    -> PNG rasters of the figures (DPI=150; override with DPI=300)
#   make clean  -> remove everything in figures/
#
# Requires pdflatex with: standalone, tikz, newtxtext, booktabs, enumitem,
# caption, hyperref, geometry.  pdftoppm is needed for `make png`.

TEXDIR := tex
OUTDIR := figures
DPI    ?= 150

FIGSRC := $(wildcard $(TEXDIR)/fig-*.tex)
FIGPDF := $(patsubst $(TEXDIR)/%.tex,$(OUTDIR)/%.pdf,$(FIGSRC))
FIGPNG := $(patsubst $(OUTDIR)/%.pdf,$(OUTDIR)/%-1.png,$(FIGPDF))

DOC    := $(OUTDIR)/loopers.pdf

.PHONY: all figs png clean
.SECONDARY: $(FIGPDF)

all: $(DOC)

figs: $(FIGPDF)

png: $(FIGPNG)

# Standalone TikZ figures. Built inside $(OUTDIR) so the aux files land there.
$(OUTDIR)/%.pdf: $(TEXDIR)/%.tex | $(OUTDIR)
	pdflatex -interaction=nonstopmode -halt-on-error \
	         -output-directory=$(OUTDIR) $<
	@rm -f $(OUTDIR)/$*.aux $(OUTDIR)/$*.log $(OUTDIR)/$*.out

# The note. Two passes for the table of contents and cross-references.
$(DOC): $(TEXDIR)/loopers.tex $(FIGPDF) | $(OUTDIR)
	pdflatex -interaction=nonstopmode -halt-on-error \
	         -output-directory=$(OUTDIR) $<
	pdflatex -interaction=nonstopmode -halt-on-error \
	         -output-directory=$(OUTDIR) $<
	@rm -f $(OUTDIR)/loopers.aux $(OUTDIR)/loopers.log \
	       $(OUTDIR)/loopers.out $(OUTDIR)/loopers.toc

$(OUTDIR)/%-1.png: $(OUTDIR)/%.pdf
	pdftoppm -png -r $(DPI) $< $(OUTDIR)/$*

$(OUTDIR):
	mkdir -p $(OUTDIR)

clean:
	rm -rf $(OUTDIR)
