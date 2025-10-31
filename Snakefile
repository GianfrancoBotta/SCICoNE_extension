from snakemake.utils import min_version
import os

configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

min_version(config["snakemake_min_version"])

container: f"docker://condaforge/mambaforge:{config['mambaforge_version']}" # to allow reproducibility

include: "workflow/rules/prepare.smk"
include: "workflow/rules/mappability.smk"
include: "workflow/rules/gc_bias.smk"
include: "workflow/rules/depth.smk"
include: "workflow/rules/qc.smk"

rule all:
    input:
        "/cluster/work/bewi/members/gbotta/exploratory_da/plots/qc_sc"