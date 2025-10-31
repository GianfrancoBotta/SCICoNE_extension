configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

# Get depth per cell in each bin
rule get_binned_coverage_and_depth_per_cell:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed",
        cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    output:
        bed_binned_cov_depth_cell=temp("/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth_sc/{sample}/{chr_region}.binned.cov_depth.bed"),
        bed_binned_cov_depth_cell_wide="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth_sc/{sample}/{chr_region}.binned.cov_depth.wide.bed"
    params:
        fraction_overlap=config["fraction_overlap"]
    log:
        "logs/qc/{chr_region}/coverage_depth/{sample}_sc.log"
    conda:
        "../envs/standalone/bedtools_sc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.bed_binned_cov_depth_cell}) && touch {output.bed_binned_cov_depth_cell};
        while read CELL; do
            samtools view -b -h -d RG:$CELL {input.bam} | \
            bedtools coverage -F {params.fraction_overlap} -a {input.bed_binned} -b - -header | \
            awk -v barcode="$CELL" 'BEGIN{{OFS="\t"}} {{bin_id=$1 ":" $2 "-" $3; print $0, barcode, bin_id}}' | cut -f 4,8,9  >> {output.bed_binned_cov_depth_cell} 2>> {log};
        done < {input.cb_file}
        mlr --tsv --implicit-csv-header reshape -s 3,1 {output.bed_binned_cov_depth_cell} > {output.bed_binned_cov_depth_cell_wide} 2>> {log}
        """

rule correct_depth_sc:
    input:
        cov_depth_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth_sc/{sample}/{chr_region}.binned.cov_depth.wide.bed",
        gc_bias_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/{chr_region}.binned.gc_bias.wide.bed",
        gc_pct_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_pct/{chr_region}.binned.gc_pct.bed"
    output:
        norm_depth_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.wide.bed"
    params:
        corr_vars=["bin_length", "gc_bias"]
    log:
        "logs/correct/{chr_region}/{sample}_correct.log"
    conda:
        "../envs/py/qc.yaml"
    script:
        "../scripts/py/correct_sc.py"

rule aggregate_depth_per_sample:
    input:
        norm_depth_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.wide.bed"
    output:
        aggr_bin_norm_depth_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.aggr_bin.bed",
        aggr_cell_norm_depth_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.aggr_cell.bed",
    log:
        "logs/correct/{chr_region}/{sample}_aggregate.log"
    conda:
        "../envs/py/qc.yaml"
    script:
        "../scripts/py/aggregate_depth_sc.py"