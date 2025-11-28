import pandas as pd

config = {
    "del_variants": "/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/sv_subset/gnomad.v4.1.sv.pass.DEL.tsv",
    'shet_file': '/hpc/users/seress01/projects/svmurt/data/s_het/shet.tsv',
    'output_dir': "/hpc/users/seress01/projects/svmurt/output/analysis/"
}

#columns are #CHROM    POS                              ID QUAL        AF SVTYPE  CHR2    END PREDICTED_LOF
df = pd.read_csv(config["del_variants"], sep="\t", names=['chrom', 'pos', 'id', 'qual', 'af', 'svtype', 'chr2', 'end', 'predicted_lof'], header=0)
#remove the rows with missing predicted_lof == .
df = df[df['predicted_lof'] != '.']

shet_df = pd.read_csv(config['shet_file'], sep='\t', names=['gene', 'hgnc', 'chrom', 'obs_lof', 'exp_lof', 'prior_mean', 'post_mean', 'post_lower_95', 'post_upper_95'])
shet_df = shet_df[['gene', 'post_mean']].rename(columns={'post_mean': 's_het'})

#make a list of all the genes in the predicted_lof column by splitting and write in a one column file without header and duplicates
lof_string = ",".join(df["predicted_lof"].astype(str).str.strip())
gene_list = sorted(set(lof_string.split(",")))
with open(f'{config["output_dir"]}lof_genes.txt', 'w') as f:
    for gene in gene_list:
        f.write(f'{gene}\n')

#python /hpc/users/seress01/tools/gene_name_mapping/gene_map.py --GNs /hpc/users/seress01/projects/svmurt/output/analysis/lof_genes.txt --output /hpc/users/seress01/projects/svmurt/output/analysis/maped_lof_genes.txt --alias_file /sc/arion/projects/lowtherlab/projects/cnv_calibration/alias_genes_gnomad_v4_browser.txt --mode unique

mapped_lof_df = pd.read_csv(f'{config["output_dir"]}maped_lof_genes.txt', sep='\t')
#make a dictionary from key gene and value gene_id
gene_to_id = dict(zip(mapped_lof_df['gene'], mapped_lof_df['gene_id']))

#convert predicted_lof column to gene_id using the gene_to_id dictionary
def convert_to_gene_id(predicted_lof):
    genes = str(predicted_lof).split(',')
    gene_ids = [str(gene_to_id.get(str(gene).strip(), str(gene).strip())) for gene in genes]
    return ','.join(gene_ids)

#add a column to df with the gene_ids
df['gene_ids'] = df['predicted_lof'].apply(convert_to_gene_id)
#add a column to df with the s_het values of the gene_ids
def get_shet_values(gene_ids):
    genes = gene_ids.split(',')
    shet_values = shet_df[shet_df['gene'].isin(genes)]['s_het'].dropna().tolist()
    return shet_values


df['s_het_values'] = df['gene_ids'].apply(get_shet_values)
#add a column to df that contains the highest s_het value of the s_het_values column
def get_max_shet(shet_values):
    if not shet_values:
        return None
    return max(shet_values)
df['max_s_het'] = df['s_het_values'].apply(get_max_shet)
#add a column to df that contains the mean s_het value of the s_het
def get_mean_shet(shet_values):
    if not shet_values:
        return None
    return sum(shet_values) / len(shet_values)
df['mean_s_het'] = df['s_het_values'].apply(get_mean_shet)

#make a set of all the genes in predicted_lof
all_genes = set()
for genes in df['predicted_lof']:
    gene_list = str(genes).split(',')
    all_genes.update([gene.strip() for gene in gene_list])
print(f'Total unique genes in predicted_lof: {len(all_genes)}')
#make a set of all the genes in gene_ids
all_gene_ids = set()
for gene_ids in df['gene_ids']:
    gene_id_list = str(gene_ids).split(',')
    all_gene_ids.update([gene_id.strip() for gene_id in gene_id_list])
print(f'Total unique gene_ids in gene_ids: {len(all_gene_ids)}')
#count how many genes in gene_ids has actually s_het value. Basically their gene_id was found in shet_df
genes_with_shet = set()
for gene_ids in df['gene_ids']:
    gene_id_list = str(gene_ids).split(',')
    for gene_id in gene_id_list:
        if not shet_df[shet_df['gene'] == gene_id].empty:
            genes_with_shet.add(gene_id.strip())
print(f'Total unique gene_ids with s_het value: {len(genes_with_shet)}')

#using plotnine and ggpplot style plot the distribution of max_s_het for the af < 0.0001 and 0.0001<=af<0.01 and af >= 0.01

df_plot = df.copy()
df_plot['max_s_het'] = pd.to_numeric(df_plot['max_s_het'], errors='coerce')
df_plot = df_plot[np.isfinite(df_plot['max_s_het'])]

p = (
    ggplot(df_plot, aes(x='max_s_het', color='af_category', fill='af_category'))
    # translucent normalized hist for area
    + geom_histogram(aes(y=after_stat('density')), bins=50, position='identity', alpha=0.25)
    # line over the same bins to look like a curve
    + geom_freqpoly(aes(y=after_stat('density')), bins=50, size=1)
    + theme_bw()
    + labs(title='Distribution of max s_het values by AF category',
           x='Max s_het', y='Density', fill='AF Category', color='AF Category')
    + coord_cartesian(xlim=(0, 1))
)

p.save('max___s_het_distribution_by_af_category.png', dpi=300)


###


df_plot = df.copy()
df_plot['mean_s_het'] = pd.to_numeric(df_plot['mean_s_het'], errors='coerce')
df_plot = df_plot[np.isfinite(df_plot['mean_s_het'])]

p3 = (
    ggplot(df_plot, aes(x='mean_s_het', color='af_category', fill='af_category'))
    # translucent normalized hist for area
    + geom_histogram(aes(y=after_stat('density')), bins=50, position='identity', alpha=0.25)
    # line over the same bins to look like a curve
    + geom_freqpoly(aes(y=after_stat('density')), bins=50, size=1)
    + theme_bw()
    + labs(title='Distribution of mean s_het values by AF category',
           x='mean s_het', y='Density', fill='AF Category', color='AF Category')
    + coord_cartesian(xlim=(0, 1))
)

p3.save('mean___s_het_distribution_by_af_category.png', dpi=300)


###



df = pd.read_csv(config["del_variants"], sep="\t", names=['chrom', 'pos', 'id', 'qual', 'af', 'svtype', 'chr2', 'end', 'predicted_lof'], header=0)
#remove the rows with missing predicted_lof == .
df_lof = df[df['predicted_lof'] != '.']

import pandas as pd
import numpy as np
from plotnine import *

# Load
df = pd.read_csv(
    config["del_variants"],
    sep="\t",
    names=['chrom', 'pos', 'id', 'qual', 'af', 'svtype', 'chr2', 'end', 'predicted_lof'],
    header=0
)

# Filter LoF rows
df_lof = df[df['predicted_lof'] != '.']

# Combine into one tidy frame
df_all = pd.concat([
    df.assign(group='All variants'),
    df_lof.assign(group='Predicted LoF')
])

# Ensure af is numeric
df_all['af'] = pd.to_numeric(df_all['af'], errors='coerce')
df_all = df_all[np.isfinite(df_all['af']) & (df_all['af'] >= 0)]

# KDE plot
p = (
    ggplot(df_all, aes(x='af', fill='group', color='group'))
    + geom_density(alpha=0.4)
    + theme_bw()
    + labs(title='AF distribution: all deletions vs predicted LoF',
           x='Allele Frequency', y='Density',
           fill='Group', color='Group')
)

p.save('af_distribution_lof_vs_all.png', dpi=300)
