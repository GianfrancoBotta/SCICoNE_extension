import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import pyarrow
import re
import seaborn as sns

from functools import reduce
from utils import *

# Create plotdir if it does not exist
os.makedirs(snakemake.output['plotdir'], exist_ok=True)

# Get sample ids
sample_ids = [os.path.basename(os.path.dirname(path)) for path in snakemake.input['bed_binned_cov_depth']]

# Load and format dataframes
gc_pct = pl.read_csv(snakemake.input['bed_binned_gc'],
                     separator='\t')
gc_pct = gc_pct.rename({'#1_usercol':'#chr', '2_usercol':'start', '3_usercol':'end'})
gc_pct = gc_pct.rename(dict(zip(gc_pct.columns, [re.sub(r"^\d+_", "", col) for col in gc_pct.columns])))

norm_depth_df_list = [None] * len(sample_ids)
gc_pct_df_list = [None] * len(sample_ids)
for i, (path_depth, path_gc_bias, outpath) in enumerate(zip(snakemake.input['bed_binned_cov_depth'], 
                                                          [os.path.join(gc_sample, 'gc_bias_metrics.txt') for gc_sample in snakemake.input['gc_bias']], 
                                                          snakemake.output['norm_depth_binned_bed'])):
    sample_id = sample_ids[i]
    # Load and format dataframes
    cov_depth = pl.read_csv(path_depth,
                    separator='\t',
                    has_header=False)
    cov_depth.columns = gc_pct.columns[0:3] + ['depth', 'covered_bases', 'bin_length', 'coverage']
    gc_bias = pl.read_csv(
        path_gc_bias,
        separator="\t",
        skip_rows=6,
        n_rows=101
    )
    gc_bias = (
        gc_bias
        .select(['GC', 'NORMALIZED_COVERAGE'])
        .rename({'GC': 'pct_gc', 'NORMALIZED_COVERAGE': 'gc_bias'})
    )
    
    # Filter out all regions with 0 coverage in off-target reads
    if snakemake.params['off_target']:
        cov_depth_filtered = cov_depth.filter(pl.col('depth') != 0)
    else:
        cov_depth_filtered = cov_depth
        
    # Join data
    cov_depth_filtered = (
        cov_depth_filtered
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
    
    gc_pct_filtered = (
        gc_pct
        .select(['#chr', 'start', 'end', 'pct_gc'])
        .join(
            cov_depth_filtered.select(['#chr', 'start', 'end']),
            on=['#chr', 'start', 'end'],
            how="inner"
            )
        )
        
    # Correct depth
    cov_depth_filtered = correct_depth(cov_depth_filtered,
                                       corr_vars=snakemake.params['corr_vars'],
                                       outpath=outpath)
    norm_depth_name = '_'.join(snakemake.params['corr_vars'] + ['norm_depth'])
    
    # Add dataframes to lists
    norm_depth_df_list[i] = cov_depth_filtered.select(['#chr', 'start', 'end', norm_depth_name]).rename({norm_depth_name: sample_id})
    gc_pct_df_list[i] = gc_pct_filtered.rename({'pct_gc': sample_id})

# Collapse list to one unique dataframe
if snakemake.params['off_target']:
    # Find common non-zero bins for all of the samples
    common_coords = reduce(
        lambda a, b: a.join(b, on=["#chr", "start", "end"], how="inner"),
        [df.select(["#chr", "start", "end"]).unique() for df in norm_depth_df_list]
    )
    # Select the common bins
    norm_depth_df_filtered_list = [
        df.join(common_coords, on=["#chr", "start", "end"], how="inner")
        .sort(["#chr", "start", "end"])
        for df in norm_depth_df_list
    ]
    gc_pct_df_filtered_list = [
        df.join(common_coords, on=["#chr", "start", "end"], how="inner")
        .sort(["#chr", "start", "end"])
        for df in gc_pct_df_list
    ]
else:
    norm_depth_df_filtered_list = norm_depth_df_list
    gc_pct_df_filtered_list = gc_pct_df_list
    
norm_depth_df = norm_depth_df_filtered_list[0]
for df in norm_depth_df_filtered_list[1:]:
    norm_depth_df = norm_depth_df.join(df, on=['#chr', 'start', 'end'], how='left')
norm_depth_df_long = norm_depth_df.unpivot(
    index=["#chr", "start", "end"],
    on=sample_ids,
    variable_name="sample_id",
    value_name="norm_depth"
)
gc_pct_df = gc_pct_df_filtered_list[0]
for df in gc_pct_df_filtered_list[1:]:
    gc_pct_df = gc_pct_df.join(df, on=['#chr', 'start', 'end'], how='left')
gc_pct_df_long = gc_pct_df.unpivot(
    index=["#chr", "start", "end"],
    on=sample_ids,
    variable_name="sample_id",
    value_name="pct_gc"
)

############ PLOTS ############
### GC pct vs normalized depth per bin
# Regress and plot GC percentage with normalized depth per bin
model = regress(df1=gc_pct_df_long,
                var1='pct_gc', 
                df2=norm_depth_df_long, 
                var2='norm_depth',
                xtransform=None,
                ytransform=np.log1p,
                xlabel='GC percentage per bin',
                ylabel='Log1p-normalized depth per bin',
                legend_loc='upper center',
                plotdir=snakemake.output['plotdir'])

### Find outliers, plot biotype, and save them in a bed file
# Load and clean gtf file
gtf_df = gtf2df(snakemake.input['gtf'])
# Find outliers using IQR
Q1 = np.quantile(norm_depth_df_long['norm_depth'], 0.25)
Q3 = np.quantile(norm_depth_df_long['norm_depth'], 0.75)
IQR = Q3 - Q1
upper_bound = Q3 + 3 * IQR

norm_depth_outliers = norm_depth_df_long.filter(np.log1p(pl.col('norm_depth')) > upper_bound)
outliers_list = (norm_depth_outliers['start'] + 1).cast(str).to_list()

# Filter gtf_df for matching outliers' regions
gtf_filtered_df = (
    gtf_df
    .select(snakemake.params['gtf_extract_vars'])
    .filter(pl.col('start').cast(pl.Utf8).is_in(outliers_list))
    .unique(subset=['seqname', 'start', 'end'])
)

# Save outliers' regions to tsv
os.makedirs(os.path.dirname(snakemake.output['outliers_tsv']), exist_ok=True)
norm_depth_outliers.write_csv(snakemake.output['outliers_tsv'], separator='\t')

# Count entries for each gene_type
for var in snakemake.params['gtf_plot_vars']:
    counts = (
        gtf_filtered_df
        .group_by(var)
        .agg(pl.len().alias("count"))
        .sort("count", descending=True)
    )
    
    # Plot barplot
    plt.figure(figsize=(10, 6))
    sns.barplot(data=counts, x=var, y="count", width=0.5)
    plt.xticks(rotation=45, ha='right')
    plt.xlabel(var)
    plt.ylabel("Frequency")
    plt.tight_layout()
    plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], f'outliers_count_{var}.png'),
        transparent=False,
        dpi=300,
        format='png'
    )
    plt.close()

### Mean log-normalized depth per sample
# Compute CV of normalized depth for each sample
depth = np.array(norm_depth_df.select(sample_ids))
depth_std = np.std(depth, axis=0)
depth_mean = np.mean(depth, axis=0)

# Plot mean normalized depth per sample and standard deviation as error bar
plt.figure(figsize=(8, 5))
plt.bar(sample_ids, depth_mean, yerr=depth_std, capsize=4)
plt.xlabel("Sample ID")
plt.ylabel("Mean normalized depth")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(
            fname=os.path.join(snakemake.output['plotdir'], 'mean_std_norm_depth_per_sample.png'),
            transparent=False,
            dpi=300,
            format='png'
        )
plt.close()

### Log-normalized depth per bin in all samples
# Plot log-normalized depth per bin
depth_df = norm_depth_df.select(sample_ids).to_pandas()
plt.figure(figsize=(16, 5))
for i, sample_id in enumerate(sample_ids):
    plt.plot(
        np.arange(depth_df.shape[0]),
        np.log1p(depth_df.iloc[:, i]),
        label=sample_id, 
        linewidth=1.0,
        alpha=0.8
    )
plt.xlabel("Bin")
plt.ylabel("Log1p-normalized depth")
plt.legend(frameon=False, bbox_to_anchor=(0.95, 1), loc='upper left')
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'log_norm_depth_per_bin.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()

### Normalized depth correlation across samples
# Compute correlations across samples' normalized depth
sp_corr_df = pd.DataFrame(np.corrcoef(depth, rowvar=False))

# Plot heatmap with correlations across samples' normalized depth
plt.figure()
sns.heatmap(sp_corr_df, cmap="vlag", center=0, annot=True, fmt=".2f")
plt.title("Normalized depth's spearman correlation across samples")
plt.xlabel("Sample")
plt.ylabel("Sample")
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'norm_depth_correlations_heatmap.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()

### Log-normalized depth per chromosome in all samples
# Compute average log-normalized depth per chromosome
chr_mean_norm_depth = (
    norm_depth_df
    .group_by("#chr")
    .agg([pl.col(s).mean().alias(s) for s in sample_ids])
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

# Plot stacked bar plot for all samples
chromosomes = chr_mean_norm_depth["#chr"].to_list()
data = np.log1p(np.array([chr_mean_norm_depth[s].to_numpy() for s in sample_ids]))
colors = plt.cm.tab10.colors
plt.figure(figsize=(12, 5))
bottom = np.zeros(len(chromosomes))
for i, sample in enumerate(sample_ids):
    plt.bar(
        chromosomes,
        data[i],
        bottom=bottom,
        label=sample,
        color=colors[i % len(colors)]
    )
    bottom += data[i]
    
plt.xlabel("Chromosome")
plt.ylabel("Mean log1p-normalized depth")
plt.xticks(rotation=45, ha="right")
plt.legend(frameon=False, bbox_to_anchor=(0.95, 1), loc='upper right')
plt.tight_layout()
plt.savefig(
        fname=os.path.join(snakemake.output['plotdir'], 'mean_log_norm_depth_per_chr.png'),
        transparent=False,
        dpi=300,
        format='png'
        )
plt.close()