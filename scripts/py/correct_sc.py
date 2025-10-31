import numpy as np
import polars as pl
import re

from utils import *

# Load data and bring it to long format
depth_sc = pl.read_csv(
    snakemake.input['cov_depth_bed'],
    separator='\t'
)
depth_sc = depth_sc.rename({'2': 'barcode'})

depth_sc_long = (
    depth_sc
    .unpivot(
        index="barcode",
        variable_name="bin",
        value_name="depth"
    )
    .with_columns(
        pl.col("bin")
        .str.split_exact(":", 1)
        .struct.rename_fields(["#chr", "bin_edges"])
        .alias("fields")
    ).unnest("fields")
    .with_columns(
        pl.col("bin_edges")
        .str.split_exact("-", 1)
        .struct.rename_fields(["start", "end"])
        .alias("fields")
    ).unnest("fields")
    .drop('bin_edges')
    .with_columns(
        pl.col("start").cast(pl.Int64).alias("start"),
        pl.col("end").cast(pl.Int64).alias("end")
    )
    .with_columns(
        (pl.col("end") - pl.col("start")).alias('bin_length')
    )
)

gc_pct = pl.read_csv(snakemake.input['gc_pct_bed'],
                     separator='\t')
gc_pct = gc_pct.rename({'#1_usercol':'#chr', '2_usercol':'start', '3_usercol':'end'})
gc_pct = gc_pct.rename(dict(zip(gc_pct.columns, [re.sub(r"^\d+_", "", col) for col in gc_pct.columns])))

gc_bias_sc = pl.read_csv(
    snakemake.input['gc_bias_bed'],
    separator='\t'
)

gc_bias_sc_long = (
    gc_bias_sc
    .unpivot(
        index="pct_gc",
        variable_name="barcode",
        value_name="gc_bias"
    )
)

# Join data
depth_sc_long_full = (
    depth_sc_long
    .join(
        gc_pct.select(['#chr', 'start', 'end', 'pct_gc']),
        on=['#chr', 'start', 'end'],
        how="left"
        )
    .with_columns(
        (pl.col("pct_gc").round(2) * 100).cast(pl.Int64).alias("pct_gc")
        )
    .join(
        gc_bias_sc_long.select(['pct_gc', 'barcode', 'gc_bias']),
        on=['pct_gc', 'barcode'],
        how='left'
        )
    )

# Correct depth
cov_depth = correct_depth(depth_sc_long_full,
                          corr_vars=snakemake.params['corr_vars'])

depth_sc_wide_full = (
    depth_sc_long_full
    .pivot(
        on='bin',
        index='barcode',
        values='depth'
    )
)

# Save wide format to tsv
os.makedirs(os.path.dirname(snakemake.output['norm_depth_bed']), exist_ok=True)
depth_sc_wide_full.write_csv(snakemake.output['norm_depth_bed'], separator='\t')

