import polars as pl
import re
import seaborn as sns

from utils import *
from scipy.stats import gaussian_kde

# Load data
gc_pct = pl.read_csv(snakemake.input['bed_binned_gc'],
                     separator='\t')
gc_pct = gc_pct.rename({'#1_usercol':'#chr', '2_usercol':'start', '3_usercol':'end'})
gc_pct = gc_pct.rename(dict(zip(gc_pct.columns, [re.sub(r"^\d+_", "", col) for col in gc_pct.columns])))

depth_sc = pl.read_csv(snakemake.input['bed_binned_cov_depth_cell_wide'],
                       separator='\t'
)
depth_sc = depth_sc.rename({'2': 'barcodes'})

# Filter out all regions with 0 coverage in off-target reads
barcodes = depth_sc["barcodes"]
depth_sc = depth_sc.drop("barcodes").select([
    col for col in depth_sc.drop("barcodes").columns
    if not depth_sc[col].sum() == 0
])
depth_sc = pl.concat([pl.DataFrame({"barcodes": barcodes}), depth_sc], how="horizontal")

depth_sc_long = (
    depth_sc
    .unpivot(
        index="barcodes",
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

# Join data
depth_sc_long = (
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
            gc_bias.select(['pct_gc', 'gc_bias']),
            on='pct_gc',
            how='left'
            )
        )

# Correct depth
depth_sc_long = correct_depth(depth_sc_long, corr_vars=snakemake.params['corr_vars'])

# Plot data
os.makedirs(snakemake.output['plotdir'], exist_ok=True)
depth_sc = (
    depth_sc_long
    .pivot(
        on='bin',
        index='barcodes',
        values='_'.join(snakemake.params['corr_vars'] + ['norm_depth'])
    )
)

# Filter out all regions with 0 normalized depth in off-target reads, including the GC only regions
barcodes = depth_sc["barcodes"]
depth_sc = depth_sc.drop("barcodes").select([
    col for col in depth_sc.drop("barcodes").columns
    if not depth_sc[col].sum() == 0
])
depth_sc = pl.concat([pl.DataFrame({"barcodes": barcodes}), depth_sc], how="horizontal")

# Save normalized depth to tsv
os.makedirs(os.path.dirname(snakemake.output['norm_depth_cell_wide_binned_bed']), exist_ok=True)
depth_sc.write_csv(snakemake.output['norm_depth_cell_wide_binned_bed'], separator='\t')

depth_sc = depth_sc.drop('barcodes')
cv = np.array(depth_sc.select([pl.col(c).std().alias(c) for c in depth_sc.columns]) / depth_sc.select([pl.col(c).mean().alias(c) for c in depth_sc.columns])).flatten()

plt.figure(figsize=(8, 5))
plt.hist(cv, bins=100, density=True, alpha=0.6, label='CV distribution')

kde = gaussian_kde(cv)
xs = np.linspace(cv.min(), cv.max(), 500)
plt.plot(xs, kde(xs), color='red', lw=2, label='KDE')

plt.xlabel("Coefficient of variation")
plt.ylabel("Density")
plt.legend()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], f'cv_distribution.png'),
        transparent=False,
        dpi=300,
        format='png'
    )
plt.close()

### Normalized depth correlation across cells
# Compute correlations across samples' normalized depth
sp_corr_df = pd.DataFrame(np.corrcoef(depth_sc, rowvar=True))

# Plot heatmap with correlations across samples' normalized depth
plt.figure()
sns.heatmap(sp_corr_df, cmap="vlag", center=0, annot=None)
plt.title("Normalized depth's spearman correlation across cells")
plt.tick_params(
    axis='x',
    which='both',
    bottom=False,
    top=False,
    labelbottom=False)
plt.tick_params(
    axis='y',
    which='both',
    left=False,
    right=False,
    labelleft=False)
plt.xlabel("Sample")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], f'norm_depth_correlations_heatmap.png'),
        transparent=False,
        dpi=300,
        format='png'
    )
plt.close()