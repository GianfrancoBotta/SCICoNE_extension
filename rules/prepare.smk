configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

### PREPARE FILES AND REFERENCES
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
        "../envs/standalone/samtools.yaml"
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

rule get_and_index_hg19:
    input:
        genome=config["align_genome"]
    output:
        genome="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa",
        index="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/index_files/ucsc_hg19.fa.fai"
    log:
        "logs/prepare/references/hg19_fa.log"
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        mkdir -p $(dirname {output.index});
        cp {input.genome} {output.genome} 2>> {log};
        samtools faidx {output.genome} -o {output.index} 2>> {log}
        """

# rule convert_hg19_gtf_to_bed:
#     input:
#         gtf="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/gtf_files/gencode.v19.annotation.exon_only.gtf"
#     output:
#         bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/bed_files/gencode.v19.annotation.exon_only.bed",
#     log:
#         "logs/prepare/references/hg19_bed.log"
#     conda:
#         "../envs/standalone/bedops.yaml"
#     shell:
#         """
#         convert2bed -i gtf < {input.gtf} > {output.bed} 2> {log}
#         """

rule create_on_off_target_binned_bed:
    input:
        var_sized_binned_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/binning/var_sized_bins.bed",
        on_target_bed=config["panel_bed"]
    output:
        on_target_binned_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/on_target/bed_files/on_target.binned.bed",
        off_target_binned_bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/off_target/bed_files/off_target.binned.bed"
    log:
        "logs/prepare/references/split_bed.log"
    conda:
        "../envs/standalone/bedtools.yaml"
    shell:
        """
        bedtools subtract -a {input.var_sized_binned_bed} -b {input.on_target_bed} | cut -f1-3 | sort -u > {output.off_target_binned_bed} 2>> {log};
        bedtools subtract -a {input.var_sized_binned_bed} -b {output.off_target_binned_bed} | cut -f1-3 | sort -u > {output.on_target_binned_bed} 2>> {log}
        """

# rule create_binned_beds:
#     input:
#         bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.bed"
#     output:
#         bed_binned="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed"    
#     params:
#         w_size=lambda wildcards: config["window_binning_size_target"]
#         if wildcards.chr_region == "on_target"
#         else config["window_binning_size_off_target"]
#     log:
#         "logs/prepare/references/{chr_region}/binned_bed.log"
#     conda:
#         "../envs/standalone/bedtools.yaml"
#     shell:
#         """
#         bedtools makewindows -b {input.bed} -w {params.w_size} | cut -f1-3 | sort -u > {output.bed_binned} 2> {log}
#         """

# To calculate GC bias
rule split_on_off_target_bams:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.bam",
        bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/{chr_region}/bed_files/{chr_region}.binned.bed"
    output:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/bams/{sample}/mapped.sorted.healthy_cells.{chr_region}.bam",
    threads: config["subset_threads"]
    log:
        "logs/prepare/bams/{sample}/{chr_region}.subset_bams_oot.log"
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools view -h -@ {threads} -L {input.bed} {input.bam} > {output.bam} 2> {log}
        """
