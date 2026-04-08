# GWAS (phenotype association) method registry.
#
# Derived from GEA_METHODS by filtering to methods that support phenotype
# association. Maintained as a separate import point so Phase 4's {source}
# wildcard unification has a clean seam to collapse them.

from .gea import GEA_METHODS

GWAS_METHODS = {
    name: cfg
    for name, cfg in GEA_METHODS.items()
    if cfg.get("supports_phenotypes", False)
}
