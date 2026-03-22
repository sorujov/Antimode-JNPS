# ── Anti-mode manuscript pipeline ─────────────────────────────────────────────
# Targets
#   make all       simulate (full R=500) + tables + pdf
#   make fast      simulate (quick R=50) + tables + pdf
#   make simulate  full Monte Carlo only
#   make tables    CSV → LaTeX tables only (requires results/ CSVs)
#   make pdf       compile manuscript only (requires submission/tables/)
#   make clean     remove LaTeX aux files and results CSVs

PYTHON   := python
PDFLATEX := pdflatex
TEX      := submission/antimode_GNST_manuscript.tex
TABLES   := submission/tables/tab_bias_rmse.tex \
             submission/tables/tab_rates.tex \
             submission/tables/tab_coverage.tex \
             submission/tables/tab_power.tex
CSVS     := results/table1_bias_rmse.csv \
             results/table2_rates.csv \
             results/table3_coverage.csv \
             results/table4_power.csv

.PHONY: all fast simulate simulate-array merge-chunks tables pdf clean

## Full pipeline: simulate (R=500) → tables → pdf
all: simulate tables pdf

## Quick pipeline: simulate (R=50) → tables → pdf
fast:
	$(PYTHON) src/simulate.py
	$(PYTHON) src/csv_to_tex.py
	cd submission && $(PDFLATEX) -interaction=nonstopmode antimode_GNST_manuscript.tex
	cd submission && $(PDFLATEX) -interaction=nonstopmode antimode_GNST_manuscript.tex

## Full Monte Carlo simulation (R=500, B=299) — single node
simulate:
	$(PYTHON) src/simulate.py --full

## Submit 10-chunk array job across multiple nodes (each chunk: 50 reps, 32 CPUs)
simulate-array:
	sbatch slurm_simulate_array.sh

## Merge raw replication chunks → final result tables (run after simulate-array completes)
merge-chunks:
	$(PYTHON) src/merge_chunks.py

## Regenerate LaTeX tables from CSVs (does not re-simulate)
tables: $(CSVS)
	$(PYTHON) src/csv_to_tex.py

## Compile manuscript PDF (twice for cross-references)
pdf: $(TABLES)
	cd submission && $(PDFLATEX) -interaction=nonstopmode antimode_GNST_manuscript.tex
	cd submission && $(PDFLATEX) -interaction=nonstopmode antimode_GNST_manuscript.tex

## Remove auxiliary files and simulation results (including chunk files)
clean:
	rm -f submission/*.aux submission/*.log submission/*.out submission/*.toc
	rm -f results/table*.csv results/raw_reps*.csv
