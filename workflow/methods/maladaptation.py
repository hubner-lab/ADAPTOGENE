# Maladaptation method registry.
#
# Each entry defines one maladaptation method. Adding a new method requires:
#   1. A dict entry here.
#   2. R scripts at the declared paths that read the standard input formats
#      and write the standard output formats (offset raster, site values, etc.).
#
# engine:                       dispatch key used by rule factories
# model_script:                 Rscript path relative to /pipeline/ (builds the model)
# offset_script:                Rscript path relative to /pipeline/ (computes genetic offset)
# cumimp_script:                Rscript path relative to /pipeline/ (cumulative importance plot)
# importance_script:            Rscript path relative to /pipeline/ (overall importance plot)
# supports_spatial:             method accepts spatial correction (PCNM) covariates
# supports_random_model:        method can build a neutral random model for comparison
# builds_model:                 method produces a separate model artifact (*.qs) before offset
# supports_cumulative_importance: method produces cumulative importance curves

MALADAPTATION_METHODS = {
    "gradient_forest": {
        "engine": "gradient_forest",
        "model_script": "scripts/gradient_forest_model.R",
        "offset_script": "scripts/gradient_forest_offset.R",
        "cumimp_script": "scripts/plot_gf_cumimp.R",
        "importance_script": "scripts/plot_gf_importance.R",
        "supports_spatial": True,
        "supports_random_model": True,
        "builds_model": True,
        "supports_cumulative_importance": True,
    },
    "geometric_offset": {
        "engine": "geometric_offset",
        "model_script": None,
        "offset_script": "scripts/geometric_offset.R",
        "cumimp_script": None,
        "importance_script": None,
        "supports_spatial": False,          # genetic.gap has no spatial argument
        "supports_random_model": False,
        "builds_model": False,              # one-call: model + offset in single script
        "supports_cumulative_importance": False,
    },
    "rda_offset": {
        "engine": "rda_offset",
        "model_script": None,
        "offset_script": "scripts/rda_offset.R",
        "cumimp_script": None,
        "importance_script": None,          # single-call script writes overall_importance.png itself
        "supports_spatial": False,          # B7 (docs/rda_research.md) — offset RDA drops Condition() by default
        "supports_random_model": False,
        "builds_model": False,              # one-call: second RDA fit + offset in single script
        "supports_cumulative_importance": False,
    },
    # Same engine and script as rda_offset, separate registry entry so BOTH the
    # uncorrected (condition_pcs: 0, canonical per B7) and the Condition()-corrected
    # fit can be produced in ONE Snakemake invocation. Previously the corrected pass
    # needed a second full invocation over every garden with its own config file,
    # which measured 15-25% of the journal-09 sweep.
    "rda_offset_corrected": {
        "engine": "rda_offset",
        "model_script": None,
        "offset_script": "scripts/rda_offset.R",
        "cumimp_script": None,
        "importance_script": None,
        "supports_spatial": False,
        "supports_random_model": False,
        "builds_model": False,
        "supports_cumulative_importance": False,
    },
}
