# ADAPTOGENE Pipeline
# vim: filetype=python
#
# Modular Snakemake workflow for population genomics.
# Config, paths, and helpers are in workflow/rules/common.smk.
# Rules are organized by pipeline mode in separate .smk files.

include: "workflow/rules/common.smk"
include: "workflow/rules/processing.smk"
include: "workflow/rules/structure.smk"
include: "workflow/rules/structure_k.smk"
include: "workflow/rules/association.smk"
if PHENO_ASSOC_CONFIGS:
    include: "workflow/rules/phenotype_assoc.smk"
if PHENO_ASSOC_CONFIGS and ASSOC_CONFIGS:
    include: "workflow/rules/overlapping.smk"
include: "workflow/rules/regionplot.smk"
include: "workflow/rules/maladaptation.smk"
include: "workflow/rules/haplotype.smk"
include: "workflow/rules/summary.smk"

rule all:
    input: get_targets(MODE)
