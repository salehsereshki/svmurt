import pandas as pd

configs = {
    'gencode_annotation': '/sc/arion/projects/lowtherlab/data/assembly/gencode.v49.primary_assembly.basic.annotation.gtf',
    'annot_types': '/hpc/users/seress01/projects/svmurt/data/input/annotation_types',
    'output_folder': '/hpc/users/seress01/projects/svmurt/data/analysis/'
}

def load_gencode_annotation(gtf_address):
    gtf_df = pd.read_csv(gtf_address, sep='\t', comment='#', header=None,
                         names=['chr', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    return gtf_df


df = load_gencode_annotation(configs['gencode_annotation'])
df = df[df['feature'] == 'gene'].copy()
df['gene_type'] = df['attribute'].apply(lambda x: [item for item in x.split('; ') if item.startswith('gene_type')][0].split(' ')[1].replace('"', ''))
df = df[df['gene_type'] == 'protein_coding'].copy()
#extract gene_id from attribute column
df['gene_id'] = df['attribute'].apply(lambda x: [item for item in x.split('; ') if item.startswith('gene_id')][0].split(' ')[1].replace('"', ''))
df = df[['chr', 'start', 'end', 'gene_id']]

#replace chri with i in chr column
df['chr'] = df['chr'].apply(lambda x: x.replace('chr', '') if 'chr' in x else x)

df.to_csv(f"{configs['output_folder']}gencode_protein_coding_genes.bed", sep='\t', index=False, header=False)
