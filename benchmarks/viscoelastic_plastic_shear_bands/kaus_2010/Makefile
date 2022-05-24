# Helper Makefile to generate the plots for the doc
ASPECT_BIN?=../../../build/aspect
MPI_NP?=6
INPUT_FILE?=kaus_2010_extension.prm
OUTPUT_DIR?=output_kaus_2010_extension

PYTHON?=python3

.PHONY: all
all: doc/kaus_2010_brick_result_time_evolution.svg \
 doc/kaus_2010_brick_result_strain_rate_profile.svg

doc/kaus_2010_brick_result_time_evolution.%: kaus_2010_time_evolution.py $(OUTPUT_DIR)/log.txt
	$(PYTHON) $< -o $(OUTPUT_DIR) -p $@

doc/kaus_2010_brick_result_strain_rate_profile.%: kaus_2010_analysis.py $(OUTPUT_DIR)/log.txt
	$(PYTHON) $< -o $(OUTPUT_DIR) -p $@

$(OUTPUT_DIR)/log.txt:
	mpirun -n $(MPI_NP) $(ASPECT_BIN) $(INPUT_FILE)

.PHONY: clean
clean:
	$(RM) doc/kaus_2010_brick_result_time_evolution.svg
	$(RM) doc/kaus_2010_brick_result_strain_rate_profile.svg