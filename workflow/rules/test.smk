rule runMQuadStandalone:
    '''
    rule to run MQuad to identify informative somatic mitochondrial variants.
    '''
    input:
        vcf = config["dir_standalone"] + "{sample_name}/cellSNP.cells.vcf.gz"
    output:
        passed_dp = config["dir_output"] + "standalone/mquad/{sample_name}/passed_dp.mtx",
        passed_ad = config["dir_output"] + "standalone/mquad/{sample_name}/passed_ad.mtx",
        passed_variant = config["dir_output"] + "standalone/mquad/{sample_name}/passed_variant_names.txt"
    conda: config["conda_env"]
    log:
        'logs/standalone/{sample_name}/runMQuadStandalone.log'
    params:
        # The input folder is the same as the cellsnp-lite output folder
        in_dir = config["dir_output"] + "standalone/cellsnplite/{sample_name}",
        # Define where MQuad should put its results
        out_dir = config["dir_output"] + "standalone/mquad/{sample_name}",
        p = config["set-threads"]["runMQuad"],
        minDP = config["minDP"]
    shell:
        '''
        mquad --vcfData {input.vcf} \
        -o {params.out_dir} \
        -p {params.p} \
        --minDP {params.minDP}
        '''

rule runVireoStandalone:
    '''
    Run vireoSNP across a range of expected clone numbers (N).
    '''
    input:
        passed_ad = rules.runMQuadStandalone.output.passed_ad,
        passed_dp = rules.runMQuadStandalone.output.passed_dp,
        barcodes = config["dir_standalone"] + "{sample_name}/cellSNP.samples.tsv",
        variant_names = rules.runMQuadStandalone.output.passed_variant
    output:
        clones_df = config["dir_output"] + "standalone/vireo/{sample_name}/barcodes_donor_ids.csv",
        variant_df = config["dir_output"] + "standalone/vireo/{sample_name}/variant_donor.csv",
        ELBO_df = config["dir_output"] + "standalone/vireo/{sample_name}/ELBO.csv"
    conda: config["conda_env"]
    log:
        'logs/standalone/{sample_name}/runVireo.log'
    params:
        n_donors = config["test_clones"]
    script:
        "../scripts/vireo.py"

rule testRenv:
    '''
    A simple dummy rule to test if the Renv environment activates correctly.
    '''
    output:
        "results/test/renv_test_passed.txt"
    conda: config["env_R"]
    script:
        "../test/R_code/analysis_R45/scr/snakemake/00_test_renv.R"