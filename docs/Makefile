OUT = BubbleProfiler

TEX = $(OUT).tex
PDF = $(OUT).pdf
BIB = $(OUT).bib
BBL = $(OUT).bbl

TEX_STY	:= \
	JHEP.bst

GENERATED_IMG	:= \
		figures/fubini_potential.pdf \
		figures/logarithmic_potential.pdf \
		figures/quartic_potential.pdf \
		figures/shooting.pdf \
		figures/ssm_eff_pot_T0.pdf \
		figures/ssm_eff_pot_Tc.pdf \
		figures/ssm_eff_pot_TN.pdf \
		figures/thin_wall_potential.pdf

IMG		:= \
		figures/alpha_tests_perturbative.pdf \
		figures/alpha_tests_shooting.pdf \
		figures/classes_potential.pdf \
		figures/classes_profile_guesser.pdf \
		figures/code_structure_v2.pdf \
		figures/compare.pdf \
		figures/fubini_from_interface.pdf \
		figures/logarithmic_from_interface.pdf \
		figures/nucleation_plot_lambda_m.pdf \
		figures/nucleation_plot_lambda_s.pdf \
		figures/nucleation_plot_Tc.pdf \
		figures/perturbative_eps_x.pdf \
		figures/perturbative_eps_y.pdf \
		figures/perturbative_x.pdf \
		figures/perturbative_y.pdf \
		figures/quartic_from_files.pdf \
		figures/thin_wall_from_interface.pdf \
		$(GENERATED_IMG)

TGZ_SRC		:= \
		$(TEX) \
		$(IMG) \
		$(BBL) \
		$(BIB) \
		$(TEX_STY)

TGZ		:= $(OUT).tar.gz

LATEX		:= latex
PDFLATEX 	:= pdflatex
BIBTEX		:= bibtex
PYTHON		:= python

LATEX_TMP	:= \
		$(patsubst %.pdf, %.aux, $(PDF)) \
		$(patsubst %.pdf, %.log, $(PDF)) \
		$(patsubst %.pdf, %.toc, $(PDF)) \
		$(patsubst %.pdf, %.out, $(PDF)) \
		$(patsubst %.pdf, %.spl, $(PDF)) \
		$(patsubst %.pdf, %.blg, $(PDF)) \
		$(patsubst %.pdf, %.bbl, $(PDF))

.PHONY: all arxiv clean distclean

all: $(PDF)

clean:
	-rm -f $(LATEX_TMP)

distclean: clean
	-rm -f $(PDF)

$(PDF): $(TEX) $(BIB) $(IMG) $(TEX_STY)
	$(PDFLATEX) $<
	$(BIBTEX) $(OUT)
	$(PDFLATEX) $<
	$(PDFLATEX) $<

figures/fubini_potential.pdf: figures/fubini_potential.py
	$(PYTHON) $<

figures/logarithmic_potential.pdf: figures/logarithmic_potential.py
	$(PYTHON) $<

figures/quartic_potential.pdf: figures/quartic_potential.py
	$(PYTHON) $<

figures/shooting.pdf: figures/shooting.py
	$(PYTHON) $<

figures/ssm_eff_pot_T0.pdf figures/ssm_eff_pot_Tc.pdf figures/ssm_eff_pot_TN.pdf: figures/plot_eff_potential.py
	$(PYTHON) $<

figures/thin_wall_potential.pdf: figures/thin_wall_potential.py
	$(PYTHON) $<

arxiv: $(TGZ)
	@true

$(TGZ): $(PDF)
	tar --transform 's,^,BubbleProfiler/,' -czf $@ $(TGZ_SRC)

share.zip: $(PDF)
	zip $@ $^ $(TGZ_SRC)
