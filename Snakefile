from snakemake.utils import min_version
import os

configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

min_version(config["snakemake_min_version"])

container: f"docker://condaforge/mambaforge:{config['mambaforge_version']}" # to allow reproducibility


rule all:
    input:
        # expand(
        #     "/cluster/work/bewi/members/gbotta/exploratory_da/plots/qc/{chr_region}",
        #     chr_region=config["chr_regions"]
        # ),
        # expand(
        #     "/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.single_cell.wide.bed",
        #     sample=config["sample_names"],
        #     chr_region=config["chr_regions"]
        # )
        expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/{chr_region}.binned.gc_bias.bed",
            sample=config["sample_names"],
            chr_region=config["chr_regions"]
        )
        

### PREPARE FILES AND REFERENCES FOR QC
rule get_healthy_cells:
    input:
        anno="/cluster/work/bewi/members/jgawron/bam_files/{sample}.healthy_cell_annotation.csv"
    output:
        anno="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    shell:
        r"""
        mkdir -p $(dirname {output.anno});
        awk -F, 'tolower($3)=="healthy" {{print $1}}' {input.anno} > {output.anno}
        """
        
rule subset_bam_healthy_cells:
    input:
        bam="/cluster/work/bewi/members/jgawron/bam_files/{sample}.mapped.sorted.bam",
        anno="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    output:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam"
    threads: config["subset_threads"]
    log:
        "logs/prepare/bams/{sample}/subset_bams_hc.log"
    conda:
        "envs/standalone/samtools.yaml"
    shell:
        """
        samtools view -h -b -@ {threads} -R {input.anno} {input.bam} > {output.bam} 2> {log};
        samtools index {output.bam}
        """

rule download_hg19_gtf_and_get_only_genes:
    output:
        gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.exon_only.gtf"
    params:
        gtf_link=config["hg19_gtf_link"],
        zipped_gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.gtf.gz",
        unzipped_gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.gtf"
    log:
        "logs/prepare/references/hg19_gtf.log"
    shell:
        r"""
        mkdir -p $(dirname {output.gtf});
        wget -O {params.zipped_gtf} {params.gtf_link} &> {log};
        gzip -dk {params.zipped_gtf};
        awk '$3=="exon"' {params.unzipped_gtf} > {output.gtf}
        """

rule convert_hg19_gtf_to_bed:
    input:
        gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.exon_only.gtf"
    output:
        bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/bed_files/gencode.v19.annotation.exon_only.bed",
    log:
        "logs/prepare/references/hg19_bed.log"
    conda:
        "envs/standalone/bedops.yaml"
    shell:
        """
        convert2bed -i gtf < {input.gtf} > {output.bed} 2> {log}
        """

rule create_on_off_target_bed:
    input:
        hg19_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/bed_files/gencode.v19.annotation.exon_only.bed",
        on_target_bed=config["panel_bed"]
    output:
        on_target_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/on_target/bed_files/on_target.bed",
        off_target_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/off_target/bed_files/off_target.bed"
    log:
        "logs/prepare/references/split_bed.log"
    conda:
        "envs/standalone/bedtools.yaml"
    shell:
        """
        cp {input.on_target_bed} {output.on_target_bed};
        bedtools subtract -a {input.hg19_bed} -b {input.on_target_bed} | cut -f1-3 | sort -u > {output.off_target_bed} 2> {log}
        """

rule create_binned_beds:
    input:
        bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.bed"
    output:
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed"    
    params:
        w_size=lambda wildcards: config["window_binning_size_target"]
        if wildcards.chr_region == "on_target"
        else config["window_binning_size_off_target"]
    log:
        "logs/prepare/references/{chr_region}/binned_bed.log"
    conda:
        "envs/standalone/bedtools.yaml"
    shell:
        """
        bedtools makewindows -b {input.bed} -w {params.w_size} | cut -f1-3 | sort -u > {output.bed_binned} 2> {log}
        """

rule split_on_off_target_bams:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
        bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.bed"
    output:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.{chr_region}.bam",
    threads: config["subset_threads"]
    log:
        "logs/prepare/bams/{sample}/{chr_region}.subset_bams_oot.log"
    conda:
        "envs/standalone/samtools.yaml"
    shell:
        """
        samtools view -h -@ {threads} -L {input.bed} {input.bam} > {output.bam} 2> {log}
        """

# rule get_cell_barcodes:
#     input:
#         bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
#     output:
#         cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/barcodes/{sample}/cells.txt"
#     threads: config["subset_threads"]
#     log:
#         "logs/prepare/barcodes/{sample}/cells.log"
#     conda:
#         "envs/standalone/samtools.yaml"
#     shell:
#         """
#         samtools view {input.bam} | awk 'BEGIN{{FS="\t"}} $16~/^RG:Z:/{{print $16}}' | cut -d":" -f3 | sort | uniq > {output.cb_file} 2> {log}
#         """

### OVERALL QC
# # Compute GC bias per sample
# rule get_binned_gc_bias:
#     input:
#         bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.{chr_region}.bam",
#         align_fa=config["align_genome"]
#     output:
#         dir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias/{sample}")
#     log:
#         "logs/qc/{chr_region}/gc_bias/{sample}.log"
#     conda:
#         "envs/standalone/picard.yaml"
#     shell:
#         """
#         mkdir -p {output.dir};
#         picard CollectGcBiasMetrics I={input.bam} O={output.dir}/gc_bias_metrics.txt CHART={output.dir}/gc_bias_metrics.pdf S={output.dir}/summary_metrics.txt R={input.align_fa} 2> {log}
#         """

# # Compute mappability bias per bin
# rule get_mappability_bias:
#     input:
#         genome=config["align_genome"]
#     output:
#         bedgraph="results/mappability/k150_mappability.bedGraph",
#         bigwig="results/mappability/k150_mappability.bw"
#     params:
#         k=150,
#         outdir="results/mappability/k150"
#     log:
#         "logs/qc/mappability/log"
#     conda:
#         "envs/standalone/umap.yaml"
#     threads: 8
#     shell:
#         """
#         mkdir -p {params.outdir}
#         umap --k {params.k} --seed-length {params.k} --threads {threads} \
#              --genome {input.genome} --output {params.outdir} &> {log}
#         """

# Compute GC percentage per bin
rule get_binned_gc_percentage:
    input:
        align_fa=config["align_genome"],
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed",
    output:
        bed_binned_gc="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_pct/{chr_region}.binned.gc_pct.bed",
    log:
        "logs/qc/{chr_region}/gc_pct/gc_pct.log"
    conda:
        "envs/standalone/bedtools.yaml"
    shell:
        """
        bedtools nuc -fi {input.align_fa} -bed {input.bed_binned} > {output.bed_binned_gc} 2> {log}
        """

# Compute coverage per bin
rule get_binned_coverage_and_depth:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed"
    output:
        bed_binned_cov_depth="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth/{sample}/{chr_region}.binned.cov_depth.bed"
    params:
        fraction_overlap=config["fraction_overlap"]
    log:
        "logs/qc/{chr_region}/coverage_depth/{sample}.log"
    conda:
        "envs/standalone/bedtools.yaml"
    shell:
        """
        bedtools coverage -F {params.fraction_overlap} -a {input.bed_binned} -b {input.bam} -header > {output.bed_binned_cov_depth} 2> {log}
        """

# QC overview
rule qc_control:
    input:
        bed_binned_gc="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_pct/{chr_region}.binned.gc_pct.bed",
        bed_binned_cov_depth=expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth/{sample}/{chr_region}.binned.cov_depth.bed",
            sample=config["sample_names"],
            chr_region=["{chr_region}"]
        ),
        gc_bias=expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias/{sample}",
            sample=config["sample_names"],
            chr_region=["{chr_region}"]
        ),
        gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.exon_only.gtf"
    output:
        plotdir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/plots/qc/{chr_region}"),
        norm_depth_binned_bed=expand(
            "/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.bed",
            sample=config["sample_names"],
            chr_region=["{chr_region}"]
        ),
        outliers_tsv="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/outliers/outliers.tsv"
    params:
        off_target=lambda wc: False if wc.chr_region == "on_target" else True,
        corr_vars=['gc_bias', 'bin_length'],
        gtf_extract_vars=['seqname', 'start', 'end', 'gene_id', 'strand', 'gene_type'],
        gtf_plot_vars=['seqname', 'strand', 'gene_type']
    log:
        "logs/qc/overview/{chr_region}.log"
    conda:
        "envs/py/qc.yaml"
    script:
        "scripts/py/qc_main.py"


### SINGLE-CELL QC
# Get depth per cell in each bin
rule get_binned_coverage_and_depth_per_cell:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
        bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed",
        cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/barcodes/{sample}/cells.txt"
    output:
        bed_binned_cov_depth_cell=temp("/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth/{sample}/{chr_region}.binned.cov_depth.single_cell.bed"),
        bed_binned_cov_depth_cell_wide="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth/{sample}/{chr_region}.binned.cov_depth.single_cell.wide.bed"
    params:
        fraction_overlap=config["fraction_overlap"]
    log:
        "logs/qc/{chr_region}/coverage_depth/{sample}_sc.log"
    conda:
        "envs/standalone/bedtools_sc.yaml"
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

# Get GC bias per cell
rule get_binned_gc_bias_per_cell:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.{chr_region}.bam",
        align_fa=config["align_genome"],
        cb_file="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/healthy_cells.txt"
    output:
        dir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/picard_sc")
    params:
        fraction_overlap=config["fraction_overlap"]
    log:
        "logs/qc/{chr_region}/gc_bias/{sample}_sc.log"
    conda:
        "envs/py/picard_sc.yaml"
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
        gc_bias_cell_wide_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias_sc/{sample}/{chr_region}.binned.gc_bias.bed"
    log:
        "logs/qc/{chr_region}/gc_bias/{sample}_sc_format.log"
    conda:
        "envs/py/picard_sc.yaml"
    script:
        "scripts/py/format_gc_bias_sc.py"

# QC overview per cell
rule sc_qc_control:
    input:
        bed_binned_gc="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_pct/{chr_region}.binned.gc_pct.bed",
        bed_binned_cov_depth_cell_wide="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/coverage_depth/{sample}/{chr_region}.binned.cov_depth.single_cell.wide.bed",
        gc_bias="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/{chr_region}/gc_bias/{sample}"
    output:
        plotdir=directory("/cluster/work/bewi/members/gbotta/exploratory_da/plots/qc_sc/{chr_region}/{sample}"),
        norm_depth_cell_wide_binned_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/correct/{chr_region}/{sample}/{chr_region}.binned.norm_depth.single_cell.wide.bed"
    params:
        corr_vars=['gc_bias', 'bin_length']
    log:
        "logs/qc/overview/{chr_region}_sc_{sample}.log"
    conda:
        "envs/py/qc.yaml"
    script:
        "scripts/py/qc_sc_main.py"