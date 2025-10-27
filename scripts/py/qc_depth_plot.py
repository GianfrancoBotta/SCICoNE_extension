import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

### Depth per bin (target regions)
# Load data
for i, t in enumerate(snakemake.input['tables_dir']):
    t_complete = os.path.join(t, 'target.regions.bed.gz')
    if i == 0:
        df = pd.read_csv(t_complete, sep="\t", compression="gzip", header=None)
    else:
        tmp = pd.read_csv(t_complete, sep="\t")
        col_to_add = tmp.iloc[:, snakemake.params['depth_idx']-1]
        df = pd.concat([df, col_to_add], axis=1)
        
# Get chromosome names to iterate over
chr_names = list(set(df.iloc[:,0]))

for chr in chr_names:
    # Create directories for each chromosome
    if not os.path.exists(os.path.join(snakemake.output['plot_dir'], f'{chr}')):
        os.makedirs(os.path.join(snakemake.output['plot_dir'], f'{chr}'))
    df_chr = df[(df.iloc[:,0] == chr)]
    ncols = min(5, len(snakemake.input['tables_dir']))
    nrows = math.ceil(len(snakemake.input['tables_dir']) / 5)

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4), sharey=True)
    for i, ax in enumerate(axes):
        ax.plot(np.arange(df_chr.shape[0]), df_chr.iloc[:, i+snakemake.params['depth_idx']-1])
        ax.set_title(f'{snakemake.params['samples'][i]}')
        ax.set_xlabel('Bins')
        if i == 0:
            ax.set_ylabel('Depth')
    fig.suptitle(f"Depth per bin off target - {chr}", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(
        fname=os.path.join(snakemake.output['plot_dir'], f'{chr}', 'depth_target.png'),
        transparent=False,
        dpi=300,
        format='png'
    )
    plt.close()

### Depth per bin (off-target regions)
# Load data
for i, t in enumerate(snakemake.input['tables_dir']):
    t_complete = os.path.join(t, 'off_target.regions.bed.gz')
    if i == 0:
        df = pd.read_csv(t_complete, sep="\t", compression="gzip", header=None)
    else:
        tmp = pd.read_csv(t_complete, sep="\t")
        col_to_add = tmp.iloc[:, snakemake.params['depth_idx']-1]
        df = pd.concat([df, col_to_add], axis=1)

# Plot data
for chr in chr_names:
    df_chr = df[(df.iloc[:,0] == chr)]
    ncols = min(5, len(snakemake.input['tables_dir']))
    nrows = math.ceil(len(snakemake.input['tables_dir']) / 5)

    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4), sharey=True)
    for i, ax in enumerate(axes):
        ax.plot(np.arange(df_chr.shape[0]), df_chr.iloc[:, i+snakemake.params['depth_idx']-1])
        ax.set_title(f'{snakemake.params['samples'][i]}')
        ax.set_xlabel('Bins')
        if i == 0:
            ax.set_ylabel('Depth')
    fig.suptitle(f"Depth per bin on target - {chr}", fontsize=16, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(
        fname=os.path.join(snakemake.output['plot_dir'], f'{chr}', 'depth_off_target.png'),
        transparent=False,
        dpi=300,
        format='png'
    )
    plt.close()