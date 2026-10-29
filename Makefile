SHELL := /bin/bash
R ?= Rscript
export R
export DATA_PATH ?= $(CURDIR)/data
export OUTPUT_PATH ?= $(CURDIR)/output
export REPLICATION_LIBRARY ?= $(CURDIR)/.replication/library
export APM_MAX_THREADS ?= 1
export APM_SEED ?= 9001
STATE := $(OUTPUT_PATH)/.replication/stamps
GRID := $(DATA_PATH)/processed_data/driver_panel_grid.csv
STAGES := conversion panels panel_reports evals eval_reports decomposition decomp_reports figure_1a_data_prep figure_1a
comma := ,
FORCED := $(subst $(comma), ,$(FORCE))
ifneq ($(filter-out all $(STAGES),$(FORCED)),)
$(error Unknown FORCE stage: $(filter-out all $(STAGES),$(FORCED)))
endif
.DEFAULT_GOAL := figure_1a
# Local jobs share the selected thread allocation; use Slurm for parallel branches.
.NOTPARALLEL:
.PHONY: setup FORCE_RUN vwh-all private-figure submit-vwh submit-figure submit-public-figure
setup:
	$(R) --vanilla setup.R

.PHONY: archive
archive: setup
	$(R) --vanilla tools/build_replication_archive.R

# make -n PLAN=1 lists only stages Make would execute and their dependencies.
# The Slurm submitter consumes this same graph without touching completion markers.
define stage
.PHONY: $(1)
$(1): $(STATE)/$(1).done
$(STATE)/$(1).done: $(addprefix $(STATE)/,$(addsuffix .done,$(2))) $(if $(filter all $(1),$(FORCED)),FORCE_RUN) $(3) | setup-once
ifeq ($$(PLAN),1)
	@: $(1) $(2)
else
	bash 'code/jobs/$(1).sh' '$(if $(filter all $(1),$(FORCED)),fresh,resume)'
endif
endef

.PHONY: setup-once
setup-once:
ifneq ($(PLAN),1)
	$(R) --vanilla setup.R
endif

$(eval $(call stage,conversion,,))
$(eval $(call stage,panels,conversion,))
$(eval $(call stage,panel_reports,panels,))
$(eval $(call stage,evals,panels,))
$(eval $(call stage,eval_reports,evals,))
$(eval $(call stage,decomposition,panels,))
$(eval $(call stage,decomp_reports,decomposition,))
$(eval $(call stage,figure_1a,,$(GRID)))
vwh-all: panel_reports eval_reports decomp_reports
private-figure: figure_1a
figure_1a_data_prep: $(GRID)
.PHONY: figure_1a_data_prep

# Supplied public grid is already a usable input. Only explicit private targets
# may create it; ordinary rendering never accesses riders.csv.
ifneq ($(filter private-figure figure_1a_data_prep,$(MAKECMDGOALS)),)
$(GRID): $(if $(filter all figure_1a_data_prep,$(FORCED)),FORCE_RUN) | setup-once
ifeq ($(PLAN),1)
	@: figure_1a_data_prep
else
	bash code/jobs/figure_1a_data_prep.sh fresh
endif
else
$(GRID):
	$(error Missing public driver grid: restore data/processed_data/driver_panel_grid.csv)
endif

submit-vwh:
	bash code/jobs/submit_pipeline.sh --target vwh-all $(if $(FORCE),--force '$(FORCE)')
submit-figure:
	bash code/jobs/submit_pipeline.sh --target private-figure $(if $(FORCE),--force '$(FORCE)')
submit-public-figure:
	bash code/jobs/submit_pipeline.sh --target figure_1a $(if $(FORCE),--force '$(FORCE)')
