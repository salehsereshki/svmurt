import pandas as pd
import numpy as np
from Bio import SeqIO
import os
from tqdm import tqdm
tqdm.pandas()


config = {
    'svtype':'DUP',
    'gene_annotation_add': '/sc/arion/projects/lowtherlab/data/assembly/gencode.v38.annotation.gtf',
    #'sv_file': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/sv_subset/gnomad.v4.1.sv.pass.DEL.tsv',
    'sv_file': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/sv_subset/gnomad.v4.1.sv.pass.DUP.tsv',
    'output_dir': '/hpc/users/seress01/projects/svmurt/output/analysis/',
    'output_plots': '/hpc/users/seress01/projects/svmurt/output/plots/',
}

os.makedirs(config['output_dir'], exist_ok=True)
os.makedirs(config['output_plots'], exist_ok=True)

svtp = config['svtype']

sv_df = pd.read_csv(config['sv_file'], sep='\t', header=0, names=['chrom', 'pos', 'id', 'qual', 'af', 'svtype', 'chr2', 'end', 'predicted_lof'])

annot_df = pd.read_csv(config['gene_annotation_add'], sep='\t', comment='#', header=None)
annot_df.columns = ['Chromosome', 'inf1', 'type', 'Start', 'End', 'inf2', 'strand', 'inf3', 'annot']


annot_df['gene_type'] = annot_df['annot'].str.extract(r'gene_type\s+"([^"]+)"')
annot_df = annot_df[annot_df['gene_type'] == 'protein_coding']
annot_df['gene_id'] = annot_df['annot'].str.extract(r'gene_id\s+"([^"]+)"')[0].str.replace(r'\.\d+$', '', regex=True)
exon_df = annot_df[annot_df.type == 'exon'][['Chromosome', 'Start', 'End', 'gene_id', 'type']]
#drop duplicates
exon_df = exon_df.drop_duplicates().reset_index(drop=True)

#for each row in sv_df, check all the exons in exon_df that overlap with the sv region (chrom, pos, end), make one column in sv_df with the list of (gene_id, exon_start, exon_end) that overlap with the sv.
def find_overlapping_exons(sv_row, exon_df_chr_scplitted, overlap_mode='full_exon'):
    chrom = sv_row['chrom']
    start = sv_row['pos']
    end = sv_row['end']
    exon_df = exon_df_chr_scplitted[chrom]
    if overlap_mode == 'any_of_exon':
        overlapping_exons = exon_df[~(exon_df['End'] < start) & ~(exon_df['Start'] > end)]
    elif overlap_mode == 'full_exon': #complete deletion of exon
        overlapping_exons = exon_df[(exon_df['Start'] >= start) & (exon_df['End'] <= end)]
    #overlapping_exons = exon_df[(exon_df['Chromosome'] == chrom) & (exon_df['Start'] <= end) & (exon_df['End'] >= start)]
    return [(row['gene_id'], row['Start'], row['End']) for _, row in overlapping_exons.iterrows()]

exon_df_chr_scplitted = {chrom: exon_df[exon_df['Chromosome'] == chrom].reset_index(drop=True) for chrom in exon_df['Chromosome'].unique()}

ol_mode = 'full_exon'  #'any_of_exon'  #'full_exon' --- 'full_deletions' is included in the name of saved file


sv_df['overlapping_exons'] = sv_df.progress_apply(lambda row: find_overlapping_exons(row, exon_df_chr_scplitted, overlap_mode=ol_mode), axis=1)
sv_df.to_csv(f'{config["output_dir"]}sv_exon_overlaps_{ol_mode}_{svtp}.tsv', sep='\t', index=False)
no_exon = sv_df[sv_df['overlapping_exons'].map(len) == 0].copy()
exon_overlapping_df = sv_df[sv_df['overlapping_exons'].map(len) > 0].copy()


###TEST UNIT CHECK DUPLICATES IN OVERLAPPING EXONS LIST ####
for i in range(10):
    lss = exon_overlapping_df[exon_overlapping_df['overlapping_exons'].map(len) > 10].iloc[i]['overlapping_exons']
    if len(lss) != len(set(lss)):
        print("There are duplicates in the overlapping exons list")
        print(i)
        #print the duplicates
        seen = set()
        duplicates = set()
        for item in lss:
            if item in seen:
                duplicates.add(item)
            else:
                seen.add(item)
        print(duplicates)



sv_df = pd.read_csv(f'{config["output_dir"]}sv_exon_overlaps_{ol_mode}_{svtp}.tsv', sep='\t')
sv_df['overlapping_exons'] = sv_df['overlapping_exons'].apply(lambda x: eval(x) if pd.notnull(x) else []).copy()
no_exon = sv_df[sv_df['overlapping_exons'].map(len) == 0].copy()
exon_overlapping_df = sv_df[sv_df['overlapping_exons'].map(len) > 0].copy()

sdfs = [exon_overlapping_df, exon_overlapping_df[exon_overlapping_df['af'] < 0.001], exon_overlapping_df[exon_overlapping_df['af'] >= 0.001]]
print(pd.DataFrame({'filter':['all', '<0.001', '>=0.001'],
                    'mean': [(s['end'] - s['pos']).mean() for s in sdfs],
                    'median': [(s['end'] - s['pos']).median() for s in sdfs]}).to_csv())


#exon_overlapping_df['overlapping_exons'] is a list of tuples (gene_id, exon_start, exon_end)
#make another column mapping the unique first element in exon_overlapping_df['overlapping_exons'] to the number of times they appear as the first element of the tuple in the list
from collections import Counter
exon_overlapping_df['overlapping_genes'] = exon_overlapping_df['overlapping_exons'].apply(lambda x: Counter(list([t[0] for t in x])))
exon_overlapping_df['overlapping_genes'] = exon_overlapping_df['overlapping_genes'].apply(lambda x: dict(x))
items = exon_overlapping_df['overlapping_genes'].map(lambda d: list(d.items()) if isinstance(d, dict) else [])
exon_overlapping_expanded = (exon_overlapping_df.drop(columns='overlapping_genes').assign(_items=items).explode('_items', ignore_index=True))
exon_overlapping_expanded[['overlapping_gene', 'overlap_value']] = pd.DataFrame(exon_overlapping_expanded['_items'].tolist(), index=exon_overlapping_expanded.index)
exon_overlapping_expanded = exon_overlapping_expanded.drop(columns='_items')

population_size = 63046

def plot_exon_deletion_distribution_per_gene(exon_overlapping_expanded, filtering, ol_mode, svtp):
    #make a bar plot for the number of SVs with  1, 2, 3, ... 20 and then >20 exon deletions per gene.
    counts = ([exon_overlapping_expanded[exon_overlapping_expanded['overlap_value'] == c]['id'].unique().size/population_size for c in range(1, 21)] +
              [exon_overlapping_expanded[exon_overlapping_expanded['overlap_value'] > 20]['id'].unique().size/population_size])
    labels = [*map(str, range(1, 21)), '>20']
    import matplotlib.pyplot as plt
    ax = pd.Series(counts, index=labels).plot(kind='bar', figsize=(8, 4), edgecolor='black')
    plt.subplots_adjust(wspace=0.1)
    ax.set_xlabel('Number of exon deletions per gene')
    ax.set_ylabel('Number of unique SVs per genome')
    plt.title(f'Distribution of SVs by number of exon deletions per gene (AF filter: {filtering})')
    plt.tight_layout()
    plt.savefig(f'{config["output_plots"]}sv_exon_deletion_distribution_per_gene_{ol_mode}_filter_{filtering}_{svtp}.png', dpi=300)
    plt.close()


# make a pd dataframe with first column the labels and three columns shows the counts based on different filtering on af
# filtering = 'all', '<0.001', '>0.001'
filterings = ['all', '<0.001', '>0.001']
counts_dict = {}
for filtering in filterings:
    if filtering == 'all':
        temp_df = exon_overlapping_expanded.copy()
    elif filtering == '<0.001':
        temp_df = exon_overlapping_expanded[exon_overlapping_expanded['af'] < 0.001].copy()
    elif filtering == '>0.001':
        temp_df = exon_overlapping_expanded[exon_overlapping_expanded['af'] > 0.001].copy()
    counts = ([temp_df[temp_df['overlap_value'] == c]['id'].unique().size for c in range(1, 21)] +
              [temp_df[temp_df['overlap_value'] > 20]['id'].unique().size])
    counts_dict[filtering] = counts
    plot_exon_deletion_distribution_per_gene(temp_df, filtering, ol_mode, svtp)
counts_df = pd.DataFrame(counts_dict, index=[*map(str, range(1, 21)), '>20'])




#make a list of counts in all rows in exon_overlapping_df['overlapping_genes'] as its ds is {gene_id1: count1, gene_id2: count2, ...}
all_gene_counts = []
all_gene_ids = []
for gene_count_dict in exon_overlapping_df['overlapping_genes']:
    all_gene_counts.extend(list(gene_count_dict.values()))
    all_gene_ids.extend(list(gene_count_dict.keys()))

gene_exon_dels_df = pd.DataFrame({'gene_id': all_gene_ids, 'deletion_count': all_gene_counts})

#make a bar plot for the number of unique genes with 1, 2, 3, ... 20 and then >20 exon deletions.
counts = []
for i in range(1, 22):
    if i < 21:
        counts.append(gene_exon_dels_df[gene_exon_dels_df['deletion_count'] == i]['gene_id'].unique().size)
    else:
        counts.append(gene_exon_dels_df[gene_exon_dels_df['deletion_count'] > 20]['gene_id'].unique().size)

import matplotlib.pyplot as plt
labels = list(range(1, 21)) + ['>20']
ax = pd.Series(counts, index=labels).plot(kind='bar', figsize=(8, 4), edgecolor='black')
#make the sapce between bars smaller
plt.subplots_adjust(wspace=0.1)
ax.set_xlabel('Number of exon deletions per gene')
ax.set_ylabel('Number of unique genes')
plt.title(f'Distribution of unique genes by number of exon deletions ({ol_mode})')
plt.tight_layout()
plt.savefig(f'{config["output_plots"]}gene_exon_deletion_distribution_{ol_mode}_{svtp}.png', dpi=300)
plt.close()




#Make a list with the number of svs with 1, 2, 3, ... 50 and then >50 overlapping exons
import numpy as np
import matplotlib.pyplot as plt
lengths = exon_overlapping_df['overlapping_exons'].map(len)
bins   = list(range(1, 52)) + [np.inf]           # edges: 1,2,...,51,inf  → 51 bins
labels = list(range(1, 51)) + ['>50']            # labels: 1..50, >50     → 51 labels
binned = pd.cut(lengths[lengths >= 1], bins=bins, labels=labels, right=False)
counts = pd.Series(binned).value_counts().sort_index()
counts = counts.reindex(list(range(1, 51)) + ['>50'], fill_value=0)
ax = counts.plot(kind='bar', figsize=(10, 4), edgecolor='black')
ax.set_xlabel('Number of overlapping exons')
ax.set_ylabel('Number of SVs')
plt.title(f'Distribution of SVs by number of overlapping exons ({ol_mode})')
plt.tight_layout()
plt.savefig(f'{config["output_plots"]}sv_exon_overlap_distribution_{ol_mode}_{svtp}.png', dpi=300)
plt.close()


#now for cumaltive counts:
thresholds = range(0, 51)
cum_counts = pd.Series({k: (lengths > k).sum() for k in thresholds})
cum_props = cum_counts / len(lengths)
ax = cum_counts.plot(kind='line', marker='o', figsize=(10, 4))
ax.set_xlabel('overlapping exons threshold')
ax.set_ylabel('SVs with more than threshold overlapping exons')
plt.title(f'Cumulative distribution of SVs by \n number of overlapping exons ({ol_mode})')
plt.xticks(thresholds, rotation=45, fontsize=8)
#tilt the x axis labels by 45 degrees and make them smaller
plt.grid()
plt.tight_layout()
plt.savefig(f'{config["output_plots"]}sv_exon_overlap_cumulative_distribution_{ol_mode}_{svtp}.png', dpi=300)
plt.close()



#plot a bar plot for the number of SVs with  1, 2, 3, ... 50  and then >50 overlapping exons- y axis not log scale, draw a black border around each bar




##TEST UNIT###
for ri in range(10):
    print(ri)
    sv_row = sv_df[sv_df['overlapping_exons'].map(len) > 0].iloc[ri]
    ls = sv_row['overlapping_exons']
    sv_s, sv_e, sv_chr = sv_row['pos'], sv_row['end'], sv_row['chrom']
    annot_rows = [annot_df[(annot_df.Chromosome == sv_chr) & (annot_df.type == 'exon') & (annot_df.Start == ls[i][1]) & (annot_df.End == ls[i][2])].iloc[0] for i in range(len(ls))]
    assert [annot_rows[i]['Start'] == ls[i][1] for i in range(len(ls))].count(True) == len(ls)
    assert [annot_rows[i]['End'] == ls[i][2] for i in range(len(ls))].count(True) == len(ls)
    assert [annot_rows[i]['gene_id'] == ls[i][0] for i in range(len(ls))].count(True) == len(ls)
    assert [annot_rows[i]['Start'] > sv_s for i in range(len(ls))].count(True) == len(ls)
    assert [annot_rows[i]['End'] < sv_e for i in range(len(ls))].count(True) == len(ls)

