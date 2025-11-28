import pandas as pd
from Bio import SeqIO
import numpy as np

config = {
    #'gnomad_bps': '/hpc/users/seress01/projects/bcr_lrpe/data/analysis/labeling/gnomAD/lifted_gBCR.bed',
    'gnomad_bps': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/new_v/temp/group1_AF_lt_0.01.tsv',
    'gnomad_bps': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/new_v/temp/group2_AF_0.01_to_0.05.tsv',
    'gnomad_bps': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/new_v/temp/group3_AF_gt_0.05.tsv',
    'hg38_gnome': '/sc/arion/projects/lowtherlab/data/assembly/hg38.fa',
    'AF_tag': 'gt0.05'
}


gdf = pd.read_csv(config['gnomad_bps'], sep='\t', header=0, names=['chr', 'pos', 'id', 'filter', 'af', 'chr2', 'end'])

# make two dataframes one 'chr', 'pos' and the other 'chr2', 'end'
gdf1 = gdf[['chr', 'pos']].copy()
gdf2 = gdf[['chr2', 'end']].copy().rename(columns={'chr2': 'chr', 'end': 'pos'})
# concatenate the two dataframes
gdf = pd.concat([gdf1, gdf2], ignore_index=True)
gdf = gdf[~gdf['chr'].isna()]

hg19_genome = {record.id: str(record.seq) for record in SeqIO.parse(config['hg38_gnome'], 'fasta')}
hg19_genome = {chrom: seq for chrom, seq in hg19_genome.items() if chrom in ['chr' + str(n) for n in range(1, 23)] or chrom in ['chrX', 'chrY']}

nmasked_regions = {chrom: np.zeros(len(seq), dtype=bool) for chrom, seq in hg19_genome.items()}
for chr in gdf['chr'].unique():
    nmasked_regions[chr][gdf[gdf['chr'] == chr]['pos'].values - 1] = 1

#count in bins of size 10kb number of 1s in masked regions
bin_size = 10000
gnomad_bps_in_hg19 = []
for chr, seq in hg19_genome.items():
    for pos in range(0, len(seq), bin_size):
        masked_count = np.sum(nmasked_regions[chr][pos:pos+bin_size])
        gnomad_bps_in_hg19.append((chr, pos, pos+bin_size, masked_count))


gnomad_bps_in_hg19_df = pd.DataFrame(gnomad_bps_in_hg19, columns=['chr', 'start', 'end', 'n_masked_regions'])

vcs = pd.DataFrame(gnomad_bps_in_hg19_df.n_masked_regions.value_counts())
#sort by "count" column vcs
vcs['mskd_rgn'] = vcs.index
vcs['mskd_rgn'] = vcs['mskd_rgn'].astype(int)
vcs = vcs.sort_values(by='mskd_rgn')

# Make a ggplot bar plotnine plot of mskd_rgn vs count
import plotnine as p9
plot = (p9.ggplot(vcs, p9.aes(x='mskd_rgn', y='count')) +
        p9.geom_bar(stat='identity') +
        p9.theme_minimal() +
        p9.labs(title='Masked Regions in hg19 Bins of Size 10000', x='Number of Masked Regions', y='Count') +
        p9.scale_x_continuous(breaks=range(0, vcs['mskd_rgn'].max() + 1, 10)))

print(plot)
# Save the plot
plot.save(f'masked_regions_hg19_bins_size_{bin_size}_{config["AF_tag"]}.png', dpi=300)

#vcs.to_csv('gnomad_bps_in_hg19_bins_size_5000.tsv', sep='\t', index=False)


#######################################################################



for record in vcf_reader:
    if record.FILTER not in [None, 'PASS']: continue
    info = record.INFO
    if info.get('SVTYPE') != 'DEL': continue
    af = info.get('AF')
    if af is None: continue
    chr2 = info.get('CHR2', '')
    end = info.get('END', '')
    line = f"{record.CHROM}\t{record.POS}\t{record.ID}\tPASS\t{af}\t{chr2}\t{end}\n"
    if af < 0.01:
        out1.write(line)
    elif af < 0.05:
        out2.write(line)
    else:
        out3.write(line)

for f in (out1, out2, out3): f.close()



# files 'group1_AF_lt_0.01.tsv', 'group1_AF_lt_0.01_awk.tsv', 'group2_AF_0.01_to_0.05.tsv', 'group2_AF_0.01_to_0.05_awk.tsv', 'group3_AF_gt_0.05.tsv', 'group3_AF_gt_0.05_awk.tsv'

#Make a list of file neames
file_names = [
    'group1_AF_lt_0.01_awk.tsv',
    'group2_AF_0.01_to_0.05_awk.tsv',
    'group3_AF_gt_0.05_awk.tsv'
]

dfs = [pd.read_csv(file_names[i], sep='\t') for i in range(3)]
# compare first and second dfs in each touple of dfs and check if their ID:AF columns match exactly
for i, (df1, df2) in enumerate(dfs):
    print(set(df1['ID'].unique()) == set(df2['ID'].unique()))
    mrgd = pd.merge(df1, df2, on='ID', suffixes=('_x', '_y'))
    mrgd['diff'] = mrgd['AF_x'] - mrgd['AF_y']
    print(len(mrgd[mrgd['diff'] > 0.01]))

