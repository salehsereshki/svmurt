import pandas as pd
import numpy as np
from Bio import SeqIO


config = {
    'gene_annotation_add': '/sc/arion/projects/lowtherlab/data/assembly/gencode.v38.annotation.gtf',
    'sv_file': '/sc/arion/projects/lowtherlab/data/variants/gnomad_4_sv/sv_subset/gnomad.v4.1.sv.pass.DEL.tsv',
    'hg38_gnome': '/sc/arion/projects/lowtherlab/data/assembly/hg38.fa',
    'shet_file': '/hpc/users/seress01/projects/svmurt/data/s_het/shet.tsv',
    # 'nmasked': '/hpc/users/seress01/projects/svmurt/data/GRCh38_Nmask.bed'
    'nmasked': '/hpc/users/seress01/projects/svmurt/data/PESR.encode.peri_all.repeats.delly.hg38.blacklist.sorted.bed'
    #'nmasked': '/hpc/users/seress01/projects/svmurt/data/bin_exclude.hg38.gatkcov.bed'
}

annot_df = pd.read_csv(config['gene_annotation_add'], sep='\t', comment='#', header=None)
annot_df.columns = ['Chromosome', 'inf1', 'type', 'Start', 'End', 'inf2', 'strand', 'inf3', 'annot']
annot_df['gene_type'] = annot_df['annot'].str.extract(r'gene_type\s+"([^"]+)"')
annot_df['gene_id'] = annot_df['annot'].str.extract(r'gene_id\s+"([^"]+)"')[0].str.replace(r'\.\d+$', '', regex=True)
annot_df = annot_df[(annot_df.type == 'gene') & (annot_df.gene_type == 'protein_coding')][['Chromosome', 'Start', 'End', 'strand', 'gene_id']]
print(f'{len(annot_df)} Genes in the annotation file')

shet_df = pd.read_csv(config['shet_file'], sep='\t', names=['gene', 'hgnc', 'chrom', 'obs_lof', 'exp_lof', 'prior_mean', 'post_mean', 'post_lower_95', 'post_upper_95'])
#merge annot_df and shet_df['gene', 'post_mean'] on annot_df['gene_id'] and shet_df['gene']
annot_df = pd.merge(annot_df, shet_df[['gene', 'post_mean']], how='left', left_on='gene_id', right_on='gene')
annot_df = annot_df.drop(columns=['gene'])
annot_df = annot_df.rename(columns={'post_mean': 's_het'})

# filter annot_df to only include the 10% lowest s_het genes
# annot_df = annot_df[~annot_df['s_het'].isna()]
# s_het_threshold = annot_df['s_het'].quantile(0.1)
# annot_df = annot_df[annot_df['s_het'] <= s_het_threshold].reset_index(drop=True)
# print(f'{len(annot_df)} Genes after filtering to the top 10% highest s_het genes')

hg38_genome = {record.id: str(record.seq) for record in SeqIO.parse(config['hg38_gnome'], 'fasta')}
hg38_genome = {chrom.replace('chr', ''): seq for chrom, seq in hg38_genome.items() if chrom.startswith('chr')}
hg38_genome = {'chr'+chrom: seq for chrom, seq in hg38_genome.items() if chrom.isdigit() and int(chrom) in range(1, 23) or chrom in ['X', 'Y']}

sv_df = pd.read_csv(config['sv_file'], sep='\t', header=0, names=['chrom', 'start', 'id', 'quality', 'af', 'svtype', 'chr2', 'end', 'predicted_lof'])

#sv_df = sv_df[sv_df['af'] > 0.001]
#sv_df = sv_df[(sv_df['end'] - sv_df['start']) > 1000]

nmasked = pd.read_csv(config['nmasked'], sep='\t', header=None, names=['chrom', 'start', 'end'])
nmasked = nmasked[nmasked['chrom'].isin(hg38_genome.keys())].reset_index(drop=True)

assert set(sv_df['chrom'].unique()) == hg38_genome.keys()
annot_df = annot_df[annot_df['Chromosome'].isin(hg38_genome.keys())].reset_index(drop=True)
assert set(annot_df['Chromosome'].unique()) == hg38_genome.keys()


flanking_range = 1000
bin_num = 20

#assert sv_df[sv_df['chrom'] != sv_df['chr2']].shape[0] == 0

# write an optimal function using numpy to get a locatoin (chr, pos) and make 20 intervals each of size 200, 10 before and 10 after the pos using intervals = np.linspace(A, B, 21)
def make_intervals(start, end, strand):
    pos = start if strand == '+' else end
    intervals = np.linspace(pos - flanking_range, pos + flanking_range, bin_num + 1)
    if strand == '-':
        intervals = intervals[::-1]
    return intervals

annot_df['tss_intervals'] = annot_df.apply(lambda row: make_intervals(row['Start'], row['End'], row['strand']), axis=1)
annot_df['tts_intervals'] = annot_df.apply(lambda row: make_intervals(row['End'], row['Start'], row['strand']), axis=1)

bpts_regions = {chrom: np.zeros(len(seq), dtype=bool) for chrom, seq in hg38_genome.items()}
# from sv_df, pull all the records for each chromosome and turn the bpts_regions[chrom] positions to 1
for chrom in sv_df['chrom'].unique():
    chrom_sv_df = sv_df[sv_df['chrom'] == chrom]
    bpts_regions[chrom][chrom_sv_df['start'].values - 1] = 1
    bpts_regions[chrom][chrom_sv_df['end'].values - 1] = 1


# gene_regions = {chrom: np.zeros(len(seq)) for chrom, seq in hg38_genome.items()}
# for _, row in annot_df.iterrows():
#     gene_regions[row['Chromosome']][row['Start']:row['End']] = gene_regions[row['Chromosome']][row['Start']:row['End']] + 1

annot_df['bpts_in_gene'] = annot_df.apply(lambda row: np.sum(bpts_regions[row['Chromosome']][row['Start']:row['End']]), axis=1)
nmasked['bpts_in_region'] = nmasked.apply(lambda row: np.sum(bpts_regions[row['chrom']][row['start']:row['end']]), axis=1)
for chrom in hg38_genome.keys():
    non_masked = len(hg38_genome[chrom]) - np.sum(nmasked[nmasked['chrom'] == chrom]['end'] - nmasked[nmasked['chrom'] == chrom]['start'])
    bpts = np.sum(bpts_regions[chrom])


#sum the bpts_regions in each interval in annot_df['tss_intervals'] and annot_df['tts_intervals']
for i in range(bin_num):
    annot_df[f'tss_bin_{i}'] = annot_df.apply(lambda row: np.sum(bpts_regions[row['Chromosome']][int(row['tss_intervals'][i]):int(row['tss_intervals'][i+1])]), axis=1)
    annot_df[f'tts_bin_{i}'] = annot_df.apply(lambda row: np.sum(bpts_regions[row['Chromosome']][int(row['tts_intervals'][i]):int(row['tts_intervals'][i+1])]), axis=1)

#Make a list of sum of tts and tss bins
tts_sums = [int(annot_df[f'tts_bin_{i}'].sum()) for i in range(bin_num)]
tss_sums = [int(annot_df[f'tss_bin_{i}'].sum()) for i in range(bin_num)]
