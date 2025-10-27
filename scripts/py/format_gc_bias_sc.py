import csv
import functools
import os
import polars as pl
import shutil

cells = []
with open(snakemake.input['cb_file'], 'r') as fd:
    cells = [line.strip() for line in fd]

gc_bias_df_list = [
    (pl.read_csv(
        os.path.join(snakemake.input['dir'], f'{cell}.gc_bias_metrics.txt'),
        separator="\t",
        skip_rows=6,
        n_rows=101
    )
    .select(['GC', 'NORMALIZED_COVERAGE'])
    # Rename 'NORMALIZED_COVERAGE' to the cell ID for the final wide format
    .rename({'GC': 'pct_gc', 'NORMALIZED_COVERAGE': cell})
    )
    for cell in cells
]

gc_bias_sc = functools.reduce(
    lambda left, right: left.join(right, on='pct_gc', how='left'),
    gc_bias_df_list
)

# Save wide format to tsv
os.makedirs(os.path.dirname(snakemake.output['norm_depth_cell_wide_binned_bed']), exist_ok=True)
gc_bias_sc.write_csv(snakemake.output['gc_bias_cell_wide_bed'], separator='\t')

# Remove directory with single-cell GC biases
shutil.rmtree(snakemake.input['dir'])
