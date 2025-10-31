configfile: "/cluster/work/bewi/members/gbotta/exploratory_da/config/config.yaml"

### GET MAPPABILITY OF THE GENOME

# Simulate reads
# rule simulate_illumina_reads:
#     input:
#         genome="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa"
#     output:
#         reads="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads.fq"
#     params:
#         coverage=lambda wildcards, input: (2 * config["sample_reads"] * 150 / sum(len(line.strip()) for line in open(input.genome) if not line.startswith(">")))
#     threads:
#         config["simulate_threads"]
#     log:
#         "logs/qc/complete/mappability/simulate_reads.log"
#     conda:
#         "../envs/standalone/art_modern.yaml"
#     shell:
#         """
#         art_modern \
#             --mode wgs \
#             --lc pe \
#             --i-file {input.genome} \
#             --o-fastq {output.reads} \
#             --builtin_qual_file HiSeq2500_150bp \
#             --read_len 150 \
#             --parallel {threads} \
#             --i-fcov {params.coverage} \
#             --pe_frag_dist_mean 300 \
#             --pe_frag_dist_std_dev 50 2> {log}
#         """

# rule split_illumina_reads:
#     input:
#         reads="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads.fq"
#     output:
#         reads_file1="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_1.fq.gz",
#         reads_file2="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_2.fq.gz"
#     params:
#         reads_file1="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_1.fq",
#         reads_file2="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_2.fq"
#     threads:
#         config["simulate_threads"]
#     log:
#         "logs/qc/complete/mappability/split_reads.log"
#     shell:
#         r"""
#         awk '{{if(NR%4==1) f=($0~"/1" ? "{params.reads_file1}" : "{params.reads_file2}"); print > f}}' {input.reads}
#         gzip -c {params.reads_file1} > {output.reads_file1};
#         gzip -c {params.reads_file2} > {output.reads_file2}
#         """

# Align simulated reads
# rule build_bowtie_index:
#      input:
#         genome="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/fasta_files/ucsc_hg19.fa",
#     output:
#         bt_index_folder=directory("/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/index_files")
#     threads:
#         config["align_threads"]
#     log:
#         "logs/qc/complete/mappability/build_index.log"
#     conda:
#         "../envs/standalone/bowtie.yaml"
#     shell:
#         """
#         mkdir -p {output.bt_index_folder};
#         bowtie2-build {input.genome} {output.bt_index_folder}/bowtie_index 2>> {log};
#         """

# rule align_illumina_reads:
#     input:
#         reads_file1="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_1.fq.gz",
#         reads_file2="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/simulation/sim_reads_2.fq.gz",
#         bt_index_folder="/cluster/work/bewi/members/gbotta/exploratory_da/results/prepare/references/complete/index_files"
#     output:
#         bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/alignment/sim_reads.aligned.sorted.bam"
#     threads:
#         config["align_threads"]
#     log:
#         "logs/qc/complete/mappability/align_reads.log"
#     conda:
#         "../envs/standalone/bowtie.yaml"
#     shell:
#         """
#         mkdir -p $(dirname {output.bam});
#         bowtie2 -x {input.bt_index_folder}/bowtie_index \
#                 -1 {input.reads_file1} -2 {input.reads_file2} \
#                 -p {threads} | \
#         samtools view -bS | samtools sort -o {output.bam} 2>> {log};
#         samtools index {output.bam} 2>> {log}
#         """

rule filter_uniquely_mapped_reads:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/alignment/sim_reads.aligned.sorted.bam"
    output:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/alignment/sim_reads.aligned.sorted.uniquely_mapped.bam"
    log:
        "logs/qc/complete/mappability/filter_reads.log"
    threads:
        config["subset_threads"]
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        samtools view -@ {threads} -b -h -q 30 -f 0x2 -F 0x4 -o {output.bam} {input.bam} 2>> {log};
        samtools index {output.bam}
        """

# Get different sized bins
rule get_variable_sized_bins:
    input:
        bam="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/alignment/sim_reads.aligned.sorted.uniquely_mapped.bam"
    output:
        bed="/cluster/work/bewi/members/gbotta/exploratory_da/results/qc/complete/mappalibilty/binning/var_sized_bins.bed"
    params:
        target_reads_per_bin=config["target_reads_per_bin"]
    log:
        "logs/qc/complete/mappability/binned_genome.log"
    conda:
        "../envs/py/mappability.yaml"
    script:
        "../scripts/py/get_mappability.py"