import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


### Average coverage per chromosome
# Extract sample name
table = snakemake.input['coverage_table']
sample = os.path.basename(os.path.dirname(table))
# Load data
df = pd.read_csv(table, sep="\t")

# Plot data
plt.figure()
plt.bar(df['#rname'], df['coverage'])
plt.xticks(rotation=90, ha='right')
plt.title('Coverage per chromosome')
plt.xlabel('Chromosomes')
plt.ylabel('Average coverage')
plt.tight_layout()
plt.savefig(fname=snakemake.output['coverage_plot'], 
            transparent=False,
            dpi=300,
            format='png')
plt.close()

# Remove M chromosome
df_noM = df[(df['#rname'] != 'chrM')]

# Plot data
plt.figure()
plt.bar(df_noM['#rname'], df_noM['coverage'])
plt.xticks(rotation=90, ha='right')
plt.title('Coverage per chromosome (without mitochondrial genes)')
plt.xlabel('Chromosomes')
plt.ylabel('Average coverage')
plt.tight_layout()
plt.savefig(fname=snakemake.output['coverage_plot_noM'], 
            transparent=False,
            dpi=300,
            format='png')
plt.close()