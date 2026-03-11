rule ruleUnzipBarcode:
    '''
    step needed to input the barcodes in cellsnp-lite.
    '''
    input:
        barcodes = config["dir_input"] + "{sample_name}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        barcodes = config["dir_output"] + "preprocess/{sample_name}/outs/barcodes.tsv"
    shadow:"minimal"
    shell:
        '''
        gunzip -c {input.barcodes} > {output.barcodes}
        '''

rule runCellSNPLite:
    '''
    rule to run the cellsnp-lite. perform the pileup of the reads.
    '''
    input:
        bam = config["dir_input"] + "{sample_name}/outs/possorted_genome_bam.bam",
        barcodes = rules.ruleUnzipBarcode.output.barcodes
    output:
        vcf = config["dir_output"] + "cellsnplite/{sample_name}/cellSNP.base.vcf.gz",
        ad = config["dir_output"] + "cellsnplite/{sample_name}/cellSNP.tag.AD.mtx",
        dp = config["dir_output"] + "cellsnplite/{sample_name}/cellSNP.tag.DP.mtx",
        barcodes = config["dir_output"] + "cellsnplite/{sample_name}/cellSNP.samples.tsv"
    conda: config["conda_env"]
    log:
        'logs/{sample_name}/runCellSNPLite.log'
    benchmark:
        'benchmarks/{sample_name}/runCellSNPLite.txt'
    resources:
        # mem_mb = config["set-resources"]["runCellSNPLite"],
        # cpus = config["set-threads"]["runCellSNPLite"]
    params:
        chrom = config["chrom"],
        p = config["set-threads"]["runCellSNPLite"],
        # Define the output directory path here to pass to the -O flag
        out_dir = config["dir_output"] + "cellsnplite/{sample_name}"
    shell:
        '''
        # We pass {params.out_dir} to -O, but Snakemake will wait until 
        # the files defined in {output} actually appear inside it.
        cellsnp-lite -s {input.bam} \
        -b {input.barcodes} \
        -O {params.out_dir} \
        -p {params.p} \
        --minMAF 0 \
        --minCOUNT 10 \
        --chrom {params.chrom} \
        --gzip
        '''

rule runMQuad:
    '''
    rule to run MQuad to identify informative somatic mitochondrial variants.
    '''
    input:
        vcf = rules.runCellSNPLite.output.vcf,
        ad = rules.runCellSNPLite.output.ad,
        dp = rules.runCellSNPLite.output.dp
    output:
        passed_dp = config["dir_output"] + "mquad/{sample_name}/passed_dp.mtx",
        passed_ad = config["dir_output"] + "mquad/{sample_name}/passed_ad.mtx",
        passed_variant = config["dir_output"] + "mquad/{sample_name}/passed_variant_names.txt"
    conda: config["conda_env"]
    log:
        'logs/{sample_name}/runMQuad.log'
    params:
        # The input folder is the same as the cellsnp-lite output folder
        in_dir = config["dir_output"] + "cellsnplite/{sample_name}",
        # Define where MQuad should put its results
        out_dir = config["dir_output"] + "mquad/{sample_name}",
        p = config["set-threads"]["runMQuad"],
        minDP = config["minDP"]
    shell:
        '''
        mquad -c {params.in_dir} \
        -o {params.out_dir} \
        -p {params.p} \
        --minDP {params.minDP}
        '''

rule runVireo:
    '''
    Run vireoSNP across a range of expected clone numbers (N).
    '''
    input:
        passed_ad = rules.runMQuad.output.passed_ad,
        passed_dp = rules.runMQuad.output.passed_dp,
        barcodes = rules.runCellSNPLite.output.barcodes,
        variant_names = rules.runMQuad.output.passed_variant
    output:
        clones_df = config["dir_output"] + "vireo/{sample_name}/barcodes_donor_ids.csv",
        variant_df = config["dir_output"] + "vireo/{sample_name}/variant_donor.csv",
        ELBO_df = config["dir_output"] + "vireo/{sample_name}/ELBO.csv"
    conda: config["conda_env"]
    log:
        'logs/{sample_name}/runVireo.log'
    params:
        n_donors = config["test_clones"]
    script:
        "../scripts/vireo.py"

rule plotMQuadResults:
    '''
    Generates exploratory heatmaps and boxplots summarizing the MQuad and vireoSNP outputs using R and ComplexHeatmap.
    '''
    input:
        # Vireo Master Tables
        clone_assignments = rules.runVireo.output.clones_df,
        variant_donors = rules.runVireo.output.variant_df,
        elbo_inits = rules.runVireo.output.ELBO_df,
        
        # MQuad Matrices (Need to access these directly from the MQuad output)
        passed_ad = rules.runMQuad.output.passed_ad,
        passed_dp = rules.runMQuad.output.passed_dp,
        variant_names = rules.runMQuad.output.passed_variant,
        
        # Original Barcodes
        barcodes = rules.ruleUnzipBarcode.output.barcodes
    output:
        # Define where the 5 PDF plots should be saved
        plot_elbo = config["dir_output"] + "plots/{sample_name}/01_ELBO.pdf",
        plot_ratio = config["dir_output"] + "plots/{sample_name}/01_hm_allelic_ratio.pdf",
        plot_prob = config["dir_output"] + "plots/{sample_name}/01_hm_donor_prob.pdf",
        plot_af_all = config["dir_output"] + "plots/{sample_name}/01_hm_donor_AF_all.pdf",
        plot_af_conf = config["dir_output"] + "plots/{sample_name}/01_hm_donor_AF_confident.pdf",
        plot_tree = config["dir_output"] + "plots/{sample_name}/01_lineageTree.pdf"
    conda: 
        config["env_R"] # Assuming you have a separate R environment defined
    params:
        # Pass the list of donors tested so R knows which N values to plot
        donors = config["test_clones"] 
    script:
        "../test/R_code/analysis_R45/scr/snakemake/01_plot_MQuad_output.R"