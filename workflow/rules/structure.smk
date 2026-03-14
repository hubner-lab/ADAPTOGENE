#=============================================================================
# MODULE 2: POPULATION STRUCTURE
#=============================================================================

rule snmf:
    """Run sNMF across K range."""
    input:  geno = W['geno']
    output: W['snmf']
    params: ks = K_START, ke = K_END, ploidy = PLOIDY, rep = REPEAT, mode = SNMF_PROJECT_MODE
    log:    f"{LOGDIR}structure/snmf.log"
    threads: CPU
    shell:
        """
        Rscript /pipeline/scripts/snmf.R {input.geno} \
            {params.ks} {params.ke} {params.ploidy} {params.rep} \
            {threads} {params.mode} > {log} 2>&1
        """

rule pca_plot:
    """Generate PCA and Tracy-Widom plots (also creates LEA PCA projections/eigenvalues)."""
    input:  lfmm = W['lfmm'], meta = O['metadata']
    output:
        pca = O['pca'],
        pca_svg = O['pca_svg'],
        tracy = O['tracy'],
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    params: plot_dir = f"{MOD_STRUCT}plots/", inter_dir = INTER
    log:    f"{LOGDIR}structure/pca_plot.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca.R \
            {input.lfmm} {input.meta} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule cross_entropy_plot:
    """Plot cross-entropy for K selection."""
    input:  snmf = W['snmf']
    output: O['cross_entropy']
    params: ks = K_START, ke = K_END, plot_dir = f"{MOD_STRUCT}plots/", inter_dir = INTER
    log:    f"{LOGDIR}structure/cross_entropy.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_cross_entropy.R \
            {input.snmf} {params.ks} {params.ke} \
            {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule extract_clusters:
    """Extract Q-matrix (cluster assignments) for specific K."""
    input:  snmf = W['snmf'], meta = O['metadata']
    output: clusters_table("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k
    log:    f"{LOGDIR}structure/extract_clusters_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/extract_clusters.R \
            {input.snmf} {input.meta} {params.k} {output} > {log} 2>&1
        """

rule structure_barplot:
    """Generate structure barplot for specific K."""
    input:  clusters = clusters_table("{k}")
    output: structure_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, plot_dir = f"{MOD_STRUCT}plots/", inter_dir = INTER
    log:    f"{LOGDIR}structure/structure_plot_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_structure.R \
            {input.clusters} {params.k} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule pca_structure_plot:
    """Generate PCA with structure pie charts for specific K (uses LEA PCA)."""
    input:
        clusters = clusters_table("{k}"),
        projections = W['pca_projections'],
        eigenvalues = W['pca_eigenvalues']
    output: pca_struct_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params: k = lambda wc: wc.k, plot_dir = f"{MOD_STRUCT}plots/", inter_dir = INTER
    log:    f"{LOGDIR}structure/pca_structure_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pca_structure.R \
            {input.clusters} {input.projections} {input.eigenvalues} \
            {params.k} {params.plot_dir} {params.inter_dir} > {log} 2>&1
        """

rule pop_diff_test:
    """Population differentiation test for specific K."""
    input:  snmf = W['snmf']
    output: pop_diff_plot("{k}")
    wildcard_constraints: k = r"\d+"
    params:
        k = lambda wc: wc.k,
        ploidy = PLOIDY,
        plot_dir = f"{MOD_STRUCT}plots/",
        inter_dir = INTER,
        scattermore_threshold = SCATTERMORE_THRESHOLD
    log:    f"{LOGDIR}structure/pop_diff_K{{k}}.log"
    shell:
        """
        Rscript /pipeline/scripts/plot_pop_diff.R \
            {input.snmf} {params.k} {params.ploidy} \
            {params.plot_dir} {params.inter_dir} {params.scattermore_threshold} > {log} 2>&1
        """
