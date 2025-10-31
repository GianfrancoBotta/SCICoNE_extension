import pandas as pd
import pysam

bam = pysam.AlignmentFile(snakemake.input['bam'], "rb")
target_reads_per_bin = snakemake.params['target_reads_per_bin']

bins = []
current_chr = None
current_start = None
read_counter = 0

for read in bam.fetch(until_eof=True):
    chrom = read.reference_name
    pos = read.reference_start

    # If we switch to a new chromosome
    if chrom != current_chr:
        if current_chr is not None:
            bins.append((current_chr, current_start, last_pos))
        current_chr = chrom
        current_start = pos
        read_counter = 0

    read_counter += 1
    last_pos = pos

    # Once we hit target_reads_per_bin, finalize the bin
    if read_counter >= target_reads_per_bin:
        bins.append((chrom, current_start, pos))
        current_start = pos + 1
        read_counter = 0

# Close last bin if it did not reach the target
if current_chr is not None and current_start < last_pos:
    bins.append((current_chr, current_start, last_pos))

df = pd.DataFrame(bins, columns=['chrom', 'start', 'end'])
df.to_csv(snakemake.output['bed'], sep="\t", index=False, header=False)