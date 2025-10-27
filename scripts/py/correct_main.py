# Load libraries
import polars as pl
import re
import matplotlib.pyplot as plt

from utils import *

# Load and format data
gc_bias = pl.read_csv(
    os.path.join(snakemake.input['gc_bias'], 'gc_bias_metrics.txt'),
    separator="\t",
    skip_rows=6,
    n_rows=101
)
gc_bias = (
    gc_bias
    .select(['GC', 'NORMALIZED_COVERAGE'])
    .rename({'GC': 'pct_gc', 'NORMALIZED_COVERAGE': 'gc_bias'})
)
gc_pct = pl.read_csv(snakemake.input['gc_pct'],
                     separator='\t')
gc_pct = gc_pct.rename({'#1_usercol':'#chr', '2_usercol':'start', '3_usercol':'end'})
gc_pct = gc_pct.rename(dict(zip(gc_pct.columns, [re.sub(r"^\d+_", "", col) for col in gc_pct.columns])))
cov_depth = pl.read_csv(snakemake.input['cov_depth'],
                        separator='\t',
                        has_header=False)
cov_depth.columns = gc_pct.columns[0:3] + ['depth', 'covered_bases', 'bin_length', 'coverage']

# Join data
cov_depth = cov_depth.join(
    gc_pct.select(['#chr', 'start', 'end', 'pct_gc']),
    on=['#chr', 'start', 'end'],
    how='left'
)
cov_depth = cov_depth.with_columns(
    (pl.col("pct_gc").round(2) * 100).cast(pl.Int64).alias("pct_gc")
)
cov_depth = cov_depth.join(
    gc_bias.select(['pct_gc', 'gc_bias']),
    on='pct_gc',
    how='left'
)

# Correct depth
cov_depth = correct_depth(cov_depth,
                          corr_vars=snakemake.params['corr_vars'],
                          outpath=snakemake.output['norm_depth_binned_bed'])

# Plot average normalized depth per chromosome
norm_col = '_'.join(snakemake.params['corr_vars'] + ['norm_depth'])
avg_depth = (
    cov_depth
    .group_by("#chr")
    .agg(pl.col(norm_col).mean().alias("avg_norm_depth"))
    .with_columns(
        pl.when(pl.col("#chr") == "chrX").then(pl.lit(23))
        .when(pl.col("#chr") == "chrY").then(pl.lit(24))
        .when(pl.col("#chr") == "chrM").then(pl.lit(25))
        .otherwise(
            pl.col("#chr")
            .str.replace("chr", "")
            .cast(pl.Int32, strict=False)
        )
        .alias("chr_num")
    )
    .sort("chr_num")
    .drop("chr_num")
)

# Create plotdir if it does not exist
os.makedirs(snakemake.output['plotdir'], exist_ok=True)

plt.figure()
plt.bar(avg_depth["#chr"], avg_depth["avg_norm_depth"])
plt.xlabel("Chromosome")
plt.ylabel("Average normalized depth")
plt.title("Average normalized depth per chromosome")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'avg_norm_depth_per_chr.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()








# Create aggregate plotdir if it does not exist
os.makedirs(os.path.join(snakemake.output['plotdir'], 'aggregate'), exist_ok=True)

# Plot barplot of CVs and means
plt.figure()
plt.bar(sample_ids, depth_cv_list, color='#636363')
plt.ylabel("Normalized depth's coefficient of variation (CV)")
plt.xlabel("Sample")
plt.ylim(0, max(depth_cv_list) * 1.2)
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'aggregate', 'coverage_cv_barplot.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()

plt.figure()
plt.bar(sample_ids, depth_mean_list, color='#636363')
plt.ylabel("Depth's mean")
plt.xlabel("Sample")
plt.ylim(0, max(depth_mean_list) * 1.2)
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'aggregate', 'coverage_mean_barplot.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()

# Compute correlations across samples' depth
if snakemake.params['off_target']:
    common_coords = reduce(
        lambda a, b: a.join(b, on=["#chr", "start", "end"], how="inner"),
        [df.select(["#chr", "start", "end"]).unique() for df in cov_depth_list]
    ) # Find common non-zero bins for all of the samples
    cov_depth_list_filtered = [
        df.join(common_coords, on=["#chr", "start", "end"], how="inner")
        .sort(["#chr", "start", "end"])
        for df in cov_depth_list
    ] # Select the common bins
else:
    cov_depth_list_filtered = cov_depth_list
    
depth_mtx = np.column_stack([
    df.select("norm_depth").to_numpy().ravel()
    for df in cov_depth_list_filtered
])
sp_corr_df = pd.DataFrame(np.corrcoef(depth_mtx, rowvar=False))

# Plot heatmap with correlations across samples' depth
plt.figure()
sns.heatmap(sp_corr_df, cmap="vlag", center=0, annot=True, fmt=".2f")
plt.title("Depth's spearman correlation")
plt.xlabel("Sample")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'aggregate', 'depth_correlations_heatmap.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()

### Plot depth per bin
depth_df = pd.DataFrame(depth_mtx)
depth_df.columns = sample_ids
    
ncols = 3
nrows = 1

fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4), sharey=True)
for i, ax in enumerate(axes):
    ax.plot(np.arange(depth_df.shape[0]), depth_df.iloc[:, i])
    ax.set_title(f'Normalized depth per bin')
    ax.set_xlabel('Bins')
    if i == 0:
        ax.set_ylabel('Normalized depth')
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'aggregate', 'depth_per_bin.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()    
