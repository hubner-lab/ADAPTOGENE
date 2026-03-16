# Benchmark Testing & Report Plan

## Context

Two benchmark datasets are prepared and ready:
- **Arabidopsis 1001G**: 1,135 accessions, ~60K SNPs, 5 chromosomes, FT10/FT16 phenotypes, full pipeline
- **Populus balsamifera**: 336 accessions, ~107K SNPs, GEA + maladaptation only (no phenotypes)

Goal: Run all pipeline modes on both datasets, verify correctness, compare with published ground truth, and produce a benchmark report (`docs/benchmarks.md`) with embedded figure links and placeholders for future features (GAPIT/BLINK).

## Testing Sequence

### Phase 1: Arabidopsis — Full Pipeline (10 modes)

Run inside Docker with `--configfile config_arabidopsis.yaml`:

| Step | Mode | Key Outputs to Verify | Est. Time |
|------|------|----------------------|-----------|
| 1 | `processing` | filtered VCF, normalized chr names (1-5), metadata, GFF | 2-5 min |
| 2 | `structure` | Q-matrices K2-K12, cross-entropy plot, PCA, Tracy-Widom | 30-45 min |
| 3 | `structure_K` | piemaps (4 climate vars), pop stats (Tajima's D, Pi, IBD, AMOVA), climate tables, density plots, correlation heatmap | 10-20 min |
| 4 | `association` | EMMAX+LFMM p-values for bio_1/4/12/15, Manhattan plots, sig SNPs, regions, genes | 20-40 min |
| 5 | `association_phenotypes` | EMMAX p-values for FT10/FT16, Manhattan plots, phenomaps | 15-30 min |
| 6 | `overlapping` | Combined GEA+GWAS regions, overlap pairs, Miami plot | 10-15 min |
| 7 | `regionplot` | Regional detail plots (needs regionplot config — see note) | 5-10 min |
| 8 | `maladaptation` | GF importance plots, genetic offset piemaps, future climate download | 20-40 min |
| 9 | `haplotype_scan` | Clustree plots, scan_status.tsv, selected_regions.tsv | 15-30 min |
| 10 | `haplotype` | crosshap viz, haplotype piemaps, assignment tables (needs epsilon selection after step 9) | 5-10 min |

**Note on regionplot**: Requires `regionplot.region`, `regionplot.traits`, `regionplot.method`, `regionplot.genes` in config. After step 4, pick a top region from `association/tables/regions_combined.tsv` and add to config before running.

**Note on haplotype**: After step 9, review clustree plots and set `haplotype.epsilon_selected` in config before running step 10.

### Phase 2: Populus — GEA + Maladaptation (8 modes)

Run inside Docker with `--configfile config_populus.yaml`:

| Step | Mode | Key Outputs to Verify | Est. Time |
|------|------|----------------------|-----------|
| 1 | `processing` | filtered VCF, normalized chr names, metadata | 2-5 min |
| 2 | `structure` | Q-matrices K2-K6, cross-entropy plot | 20-30 min |
| 3 | `structure_K` | piemaps (5 climate vars), pop stats, climate tables | 10-20 min |
| 4 | `association` | EMMAX+LFMM p-values for bio_1/2/4/12/15, Manhattan plots, regions, genes | 20-40 min |
| 5 | `regionplot` | Regional detail plots (needs config — same note as above) | 5-10 min |
| 6 | `maladaptation` | GF importance, genetic offset piemaps (compare with Fitzpatrick et al. common garden) | 20-40 min |
| 7 | `haplotype_scan` | Clustree plots, scan status | 15-30 min |
| 8 | `haplotype` | crosshap viz, haplotype piemaps (needs epsilon selection) | 5-10 min |

**NOT applicable for Populus**: `association_phenotypes`, `overlapping` (no phenotype data).

### Error Handling

After each mode:
1. Check exit code (non-zero = failure)
2. Check log files in `{PROJECT}_logs/{module}/`
3. Verify key output files exist and are non-empty
4. If a mode fails, diagnose from logs before proceeding

## Benchmark Report: `docs/benchmarks.md`

Structure of the report with embedded relative links to figures:

```markdown
# ADAPTOGENE Benchmark Report

## Overview
- Two datasets, N modes tested, pipeline version, date
- Summary table: dataset × mode × status (pass/fail/skip)

## Dataset 1: Arabidopsis 1001 Genomes

### Data Summary
- 1,135 accessions, 60K SNPs, 5 chromosomes
- Phenotypes: FT10 (1,003 values), FT16 (970 values)
- Source: 1001 Genomes Consortium 2016

### Population Structure
- Cross-entropy plot → optimal K
- PCA plot with population labels
- ![Cross-entropy](../Arabidopsis_results/structure/plots/cross_entropy_K2_K12.png)
- ![PCA](../Arabidopsis_results/structure/plots/pca.png)

### Climate Association (GEA)
- Manhattan plots for bio_1, bio_4, bio_12, bio_15
- Sig SNP counts per trait per method
- ![Manhattan combined](../Arabidopsis_results/association/plots/manhattan_combined_K8.png)
- Regions table (top 10): how many overlap FRI/FLC?

### Phenotype Association (GWAS) — Ground Truth Comparison
**Key validation**: Do we detect FRI (chr4:~269kb) and FLC (chr5:~3.17Mb)?

- Manhattan plots for FT10, FT16
- ![Manhattan FT10](../Arabidopsis_results/phenotype_association/plots/manhattan/EMMAX/manhattan_FT10_K8_bonf.png)
- ![Manhattan FT16](../Arabidopsis_results/phenotype_association/plots/manhattan/EMMAX/manhattan_FT16_K8_bonf.png)
- Table: sig SNPs within 500kb of FRI and FLC
- Comparison with Atwell et al. 2010 known loci

### GEA × GWAS Overlap
- Overlap summary: how many GEA regions co-locate with GWAS hits?
- Miami plot
- Biological interpretation: flowering time loci under both climate and phenotype selection

### Maladaptation (Gradient Forest)
- Variable importance: which climate predictors drive genetic offset?
- Genetic offset map: where are populations most maladapted under SSP585?
- ![GF importance](../Arabidopsis_results/maladaptation/plots/overall_importance_*.png)
- ![Genetic offset](../Arabidopsis_results/maladaptation/plots/genetic_offset_piemap_*.png)

### Haplotype Analysis
- Top regions from association: haplotype structure matches population structure?
- Clustree plot for epsilon selection

## Dataset 2: Populus balsamifera

### Data Summary
- 336 accessions, ~107K SNPs, 42 populations, North American transect

### Population Structure
- K=3 matches latitudinal gradient?
- PCA shows north-south cline?

### Climate Association (GEA)
- Manhattan plots for 5 climate variables
- Latitudinal cline SNPs detected?

### Maladaptation — Ground Truth Comparison
**Key validation**: Does genetic offset correlate with common garden fitness (Fitzpatrick et al. 2021)?

- Gradient Forest variable importance
- Genetic offset map across North American range
- Comparison with published GF results (Fitzpatrick Table 1)
- Qualitative agreement: northern populations predicted more maladapted?

### Haplotype Analysis
- GEA-derived regions: haplotype diversity along latitudinal gradient?

## Future Directions (Placeholders)

### GAPIT/BLINK Integration (Planned)
- **Status**: Not yet implemented
- **Expected improvement**: BLINK model for GWAS expected to improve power for detecting FRI/FLC in Arabidopsis FT10/FT16 due to multi-locus framework
- **Benchmark plan**: Re-run `association_phenotypes` with GAPIT/BLINK configs, compare sig SNP overlap with EMMAX results near FRI (chr4:269k) and FLC (chr5:3.17M)
- **For GEA**: GAPIT could complement LFMM for climate association — benchmark against same climate predictors

### Additional GEA Methods (Planned)
- **RONA** (Risk of Non-Adaptedness): complement Gradient Forest genetic offset
- **LEA offset**: Alternative genetic offset using latent factors
- **Benchmark plan**: Compare genetic offset predictions across methods using Populus common garden data as ground truth

### Phenotype Processing Pipeline (Planned)
- Prerequisite for GAPIT and HapMap analysis modes
- Will add phenotype QC, transformation, and BLUE/BLUP computation
- **Benchmark plan**: Compare raw vs processed phenotype GWAS results for FT10/FT16

## Methods
- Pipeline version, Docker image, Snakemake version
- Hardware specs (CPU, RAM, runtime per mode)
- Parameter settings from config YAML files
```

## Execution Plan

### Step-by-step execution commands

The executor should:

1. **Run Arabidopsis modes 1-6 sequentially** (processing through overlapping)
2. **Inspect association results** to configure regionplot (pick top region)
3. **Run regionplot** with configured region
4. **Run maladaptation**
5. **Run haplotype_scan**, inspect clustree, set epsilon
6. **Run haplotype**
7. **Repeat for Populus** (modes 1-4, then regionplot, maladaptation, haplotype_scan, haplotype)
8. **Write `docs/benchmarks.md`** by reading output tables and linking figures

### Docker command template
```bash
docker run --user $(id -u):$(id -g) --rm --memory=20g -v $PWD:/pipeline adaptogene:latest \
  snakemake -c10 -s Snakefile --config mode=MODE --configfile CONFIG.yaml --scheduler greedy
```

### Validation after each mode
```bash
# Check exit code
echo $?
# Check summary table
tail -5 {PROJECT}_results/pipeline_summary.tsv
# Check logs for errors
grep -i "error\|fail\|exception" {PROJECT}_logs/{module}/*.log
```

### Ground truth comparison (automated)
After association_phenotypes completes for Arabidopsis:
```bash
# Check for sig SNPs near FRI (chr4:250k-290k) and FLC (chr5:3.0M-3.4M)
awk -F'\t' '$1==4 && $2>250000 && $2<290000' Arabidopsis_results/phenotype_association/tables/EMMAX/EMMAX_sig_snps_bonf.tsv
awk -F'\t' '$1==5 && $2>3000000 && $2<3400000' Arabidopsis_results/phenotype_association/tables/EMMAX/EMMAX_sig_snps_bonf.tsv
```

## Files to Create/Modify

| File | Action |
|------|--------|
| `docs/benchmarks.md` | **Create** — benchmark report with embedded figures |
| `config_arabidopsis.yaml` | **Modify** — add regionplot config + haplotype epsilon after interactive steps |
| `config_populus.yaml` | **Modify** — add regionplot config + haplotype epsilon after interactive steps |

## Verification

1. All modes exit 0 for both datasets
2. `pipeline_summary.tsv` has entries for each completed mode
3. Key plots exist and are non-empty (Manhattan, PCA, piemaps, GF offset)
4. Ground truth loci detected (or documented as not detected with explanation)
5. Report renders correctly in GitHub markdown with relative image links
