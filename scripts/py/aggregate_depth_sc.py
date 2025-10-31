import os
import polars as pl

norm_depth_sc = pl.read_csv(
    snakemake.input['norm_depth_bed'],
    separator='\t'
)

aggr_per_bin = (
    norm_depth_sc
    .drop("barcode")
    .sum()
    .transpose(include_header=True)
    .rename({"column": "bin", "column_0": "norm_depth"})
)
aggr_per_cell = (
    norm_depth_sc
    .with_columns(
        pl.sum_horizontal(norm_depth_sc.drop("barcode")).alias("norm_depth")
    ).select(["barcode", "norm_depth"])
)

# Save to tsv
os.makedirs(os.path.dirname(snakemake.output['aggr_bin_norm_depth_bed']), exist_ok=True)
aggr_per_bin.write_csv(snakemake.output['aggr_bin_norm_depth_bed'], separator='\t')
aggr_per_cell.write_csv(snakemake.output['aggr_cell_norm_depth_bed'], separator='\t')