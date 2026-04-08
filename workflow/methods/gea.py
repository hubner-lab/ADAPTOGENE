# GEA (Genotype-Environment Association) method registry.
#
# Each entry defines one association method. Adding a new method requires:
#   1. A dict entry here.
#   2. An R script at the declared 'script' path that reads the standard
#      input formats for its engine and writes pvalue/qvalue TSVs.
#
# engine:              dispatch key used by rule factories ("emmax"|"lfmm"|"gapit")
# script:              Rscript path relative to /pipeline/
# gapit_model:         model name passed to gapit.R (gapit engine only)
# supports_phenotypes: method can be used in phenotype_association (GWAS) mode
# supports_drop:       method supports DROP mode (per-trait sample subsetting)

def _gapit(model):
    return {
        "engine": "gapit",
        "script": "scripts/gapit.R",
        "gapit_model": model,
        "supports_phenotypes": True,
        "supports_drop": True,
    }

GEA_METHODS = {
    "EMMAX": {
        "engine": "emmax",
        "script": "scripts/emmax.R",
        "supports_phenotypes": True,
        "supports_drop": True,
    },
    "LFMM": {
        "engine": "lfmm",
        "script": "scripts/lfmm.R",
        "supports_phenotypes": False,
        "supports_drop": False,
    },
    "BLINK":   _gapit("BLINK"),
    "FarmCPU": _gapit("FarmCPU"),
    "MLM":     _gapit("MLM"),
    "MLMM":    _gapit("MLMM"),
    "GLM":     _gapit("GLM"),
    "CMLM":    _gapit("CMLM"),
    "ECMLM":   _gapit("ECMLM"),
    "SUPER":   _gapit("SUPER"),
}
