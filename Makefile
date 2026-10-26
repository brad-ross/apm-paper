R ?= Rscript
FORCE ?=
FORCE_ARG = $(if $(strip $(FORCE)),--force $(FORCE),)
.DEFAULT_GOAL := figure-1a
.PHONY: setup setup-public figure-1a private-figure vwh-all submit-vwh submit-figure conversion panels panel-reports evals eval-reports decomposition decomp-reports
setup:
	$(R) --vanilla setup.R --scope vwh
setup-public:
	$(R) --vanilla setup.R --scope public
figure-1a private-figure vwh-all conversion panels panel-reports evals eval-reports decomposition decomp-reports:
	$(R) --vanilla run.R --target $@ --mode local $(FORCE_ARG)
submit-vwh:
	bash code/submit_pipeline.sh --target vwh-all $(FORCE_ARG)
submit-figure:
	bash code/submit_pipeline.sh --target private-figure $(FORCE_ARG)
