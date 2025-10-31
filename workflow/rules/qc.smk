configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

# QC overview
rule sc_qc_control:
    input:
        norm_depth_aggr_bin =expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{region}/{sample}/{region}.binned.norm_depth.aggr_bin.bed",
            sample=config["sample_names"],
            region=config["chr_regions"]
        ),
        norm_depth_aggr_cell =expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{region}/{sample}/{region}.binned.norm_depth.aggr_cell.bed",
            sample=config["sample_names"],
            region=config["chr_regions"]
        )
    output:
        plotdir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/plots/qc_sc")
    params:
        sample_ids=config["sample_names"],
        chr_regions=config["chr_regions"]
    log:
        "logs/qc/overview_sc.log"
    conda:
        "../envs/r/qc.yaml"
    script:
        "../scripts/r/qc_sc.R"