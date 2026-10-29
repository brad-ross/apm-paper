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
STAGES := conversion panels panel-reports evals eval-reports decomposition decomp-reports prepare-figure figure-1a
comma := ,
FORCED := $(subst $(comma), ,$(FORCE))
ifneq ($(filter-out all $(STAGES),$(FORCED)),)
$(error Unknown FORCE stage: $(filter-out all $(STAGES),$(FORCED)))
endif
.DEFAULT_GOAL := figure-1a
# Local jobs share the selected thread allocation; use Slurm for parallel branches.
.NOTPARALLEL:
.PHONY: setup setup-public setup-vwh FORCE_RUN vwh-all private-figure submit-vwh submit-figure submit-public-figure
setup: setup-vwh
setup-vwh:
	$(R) --vanilla setup.R --scope vwh
setup-public:
	$(R) --vanilla setup.R --scope public

# make -n PLAN=1 lists only stages Make would execute and their dependencies.
# The Slurm submitter consumes this same graph without touching completion markers.
define stage
.PHONY: $(1)
$(1): $(STATE)/$(1).done
$(STATE)/$(1).done: $(addprefix $(STATE)/,$(addsuffix .done,$(2))) $(if $(filter all $(1),$(FORCED)),FORCE_RUN) $(3) | setup-once
ifeq ($$(PLAN),1)
	@: $(1) $(2)
else
	bash code/jobs/run_stage.sh '$(1)' '$(if $(filter all $(1),$(FORCED)),fresh,resume)'
endif
endef

.PHONY: setup-once
setup-once:
ifneq ($(PLAN),1)
	$(R) --vanilla setup.R --scope $(if $(filter vwh-all conversion panels panel-reports evals eval-reports decomposition decomp-reports,$(MAKECMDGOALS)),vwh,public)
endif

$(eval $(call stage,conversion,,))
$(eval $(call stage,panels,conversion,))
$(eval $(call stage,panel-reports,panels,))
$(eval $(call stage,evals,panels,))
$(eval $(call stage,eval-reports,evals,))
$(eval $(call stage,decomposition,panels,))
$(eval $(call stage,decomp-reports,decomposition,))
$(eval $(call stage,figure-1a,,$(GRID)))
vwh-all: panel-reports eval-reports decomp-reports
private-figure: figure-1a
prepare-figure: $(GRID)
.PHONY: prepare-figure

# Supplied public grid is already a usable input. Only explicit private targets
# may create it; ordinary rendering never accesses riders.csv.
ifneq ($(filter private-figure prepare-figure,$(MAKECMDGOALS)),)
$(GRID): $(if $(filter all prepare-figure,$(FORCED)),FORCE_RUN) | setup-once
ifeq ($(PLAN),1)
	@: prepare-figure
else
	bash code/jobs/run_stage.sh prepare-figure fresh
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
	bash code/jobs/submit_pipeline.sh --target figure-1a $(if $(FORCE),--force '$(FORCE)')
