configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

# Compute GC percentage per bin
rule get_binned_gc_percentage:
    input:
        align_fa="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa",
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed"
    output:
        bed_binned_gc="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_pct/{chr_region}.binned.gc_pct.bed"
    log:
        "logs/qc/{chr_region}/gc_pct/gc_pct.log"
    conda:
        "../envs/standalone/bedtools.yaml"
    shell:
        """
        bedtools nuc -fi {input.align_fa} -bed {input.bed_binned} > {output.bed_binned_gc} 2> {log}
        """

# Get GC bias per cell
rule get_binned_gc_bias_per_cell:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.{chr_region}.bam",
        align_fa="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa",
        cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    output:
        dir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/picard_sc")
    params:
        fraction_overlap=config["fraction_overlap"]
    log:
        "logs/qc/{chr_region}/gc_bias/{sample}_sc.log"
    conda:
        "../envs/py/picard_sc.yaml"
    shell:
        """
        mkdir -p {output.dir};
        while read CELL; do
            samtools view -b -h -d RG:$CELL {input.bam} | \
            picard CollectGcBiasMetrics I=/dev/stdin O={output.dir}/$CELL.gc_bias_metrics.txt CHART={output.dir}/$CELL.gc_bias_metrics.pdf S={output.dir}/$CELL.summary_metrics.txt R={input.align_fa} >> {log}
        done < {input.cb_file}
        """

rule format_gc_bias_per_cell:
    input:
        dir="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/picard_sc",
        cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    output:
        bed_binned_gc_bias_cell_wide="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/{chr_region}.binned.gc_bias.wide.bed"
    log:
        "logs/qc/{chr_region}/gc_bias/{sample}_sc_format.log"
    conda:
        "../envs/py/picard_sc.yaml"
    script:
        "../scripts/py/format_gc_bias_sc.py"