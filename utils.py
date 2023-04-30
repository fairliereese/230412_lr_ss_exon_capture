import pdb
import scanpy as sc
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
import os
import numpy as np
# import pyranges as pr
# import pyranges as pyranges
import scipy.stats as st
import scipy
# import swan_vis as swan
import pdb
import pysam

def reverse_alignment(infile, outfile, threads=1):

    reverse_strand = {0: 16, 16: 0}

    if infile.endswith('.bam'):
        in_mode = 'rb'
    else:
        in_mode = 'r'
    if outfile.endswith('bam'):
        out_mode = 'wb'
    else:
        out_mode = 'w'
    input =  pysam.AlignmentFile(infile, in_mode, threads=threads)
    output = pysam.AlignmentFile(outfile, out_mode, template=input, threads=threads)

    for read in input:
        if read.has_tag('ts') and read.flag in reverse_strand:
            if read.get_tag('ts') == '-':
                read.flag = reverse_strand[read.flag]
                read.set_tag('ts', '+')

        output.write(read)
    input.close()
    output.close()

def get_biotype_map():
    """
    Get a dictionary mapping each gene type to a more general biotype
    """
    map = {'protein_coding': ['protein_coding'],
           'lncRNA': ['lincRNA',
                      'lncRNA',
                      'processed_transcript',
                      'sense_intronic',
                      '3prime_overlapping_ncRNA',
                      'bidirectional_promoter_lncRNA',
                      'sense_overlapping',
                      'non_coding',
                      'macro_lncRNA',
                      'antisense'],
           'pseudogene': ['unprocessed_pseudogene',
                          'translated_unprocessed_pseudogene',
                          'transcribed_unprocessed_pseudogene',
                          'processed_pseudogene',
                          'transcribed_processed_pseudogene',
                          'transcribed_unitary_pseudogene',
                          'unitary_pseudogene',
                          'polymorphic_pseudogene',
                          'pseudogene',
                          'translated_processed_pseudogene'],
           'miRNA': ['miRNA'],
           'other': ['snRNA', 'vault_RNA',
                     'misc_RNA', 'TEC',
                     'snoRNA', 'scaRNA',
                     'rRNA_pseudogene', 'rRNA',
                     'IG_V_pseudogene',
                     'scRNA', 'IG_V_gene',
                     'IG_C_gene', 'IG_J_gene',
                     'sRNA', 'ribozyme',
                     'vaultRNA', 'TR_C_gene',
                     'TR_J_gene', 'TR_V_gene',
                     'TR_V_pseudogene', 'TR_D_gene',
                     'IG_C_pseudogene', 'IG_D_gene',
                     'IG_pseudogene', 'Mt_tRNA',
                     'Mt_rRNA', 'TR_J_pseudogene',
                     'IG_J_pseudogene']}
    return map

def get_gene_info(gtf, o):
    df = pr.read_gtf(gtf, as_df=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # only gene
    df = df.loc[df.Feature == 'gene'].copy(deep=True)

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    # gene length
    df['length'] = (df.Start-df.End).abs()

    # add TF info
    df['tf'] = False
    d = os.path.dirname(__file__)
    tf_df = pd.read_csv('{}/../refs/biomart_tf_gids.tsv'.format(d), sep='\t')
    tf_gids = tf_df['Gene stable ID'].unique().tolist()
    df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df = df[['gid', 'gname', 'length', 'biotype', 'biotype_category', 'tf']]
    df.to_csv(o, sep='\t', index=False)

def get_transcript_info(gtf, o):
    """
    Get a file with relevant information about each transcript in a gtf

    Parameters:
        gtf (str): Path to gtf
        o (str): Output file name
    """

    df = pr.read_gtf(gtf, as_df=True)

    # remove sirvs and erccs
    print(len(df.index))
    df = df.loc[(~df.Chromosome.str.contains('SIRV'))&~(df.Chromosome.str.contains('ERCC'))]
    print(len(df.index))

    # only exons
    df = df.loc[df.Feature == 'exon'].copy(deep=True)

    # rename some columns
    m = {'gene_id': 'gid',
         'gene_name': 'gname',
         'transcript_id': 'tid',
         'gene_type': 'biotype'}
    df.rename(m, axis=1, inplace=True)

    map = get_biotype_map()

    beeps = []
    for key, item in map.items():
        beeps += item

    set(df.biotype.unique().tolist())-set(beeps)

    # pivot map
    biotype_map = {}
    for key, biotypes in map.items():
        for biotype in biotypes:
            biotype_map[biotype] = key

    # then add map to df
    df['biotype_category'] = df.biotype.map(biotype_map)

    df['exon_len'] = (df.Start-df.End).abs()+1

    df = df[['gid', 'gname', 'tid', 'exon_len', 'biotype', 'biotype_category']]
    df_copy = df[['gid', 'gname', 'tid', 'biotype', 'biotype_category']].copy(deep=True)
    df_copy = df_copy.drop_duplicates(keep='first')

    df = df.groupby('tid').sum().reset_index()
    df.rename({'exon_len': 't_len'}, axis=1, inplace=True)
    df = df.merge(df_copy, on='tid', how='left')

    # add TF info
    df['tf'] = False
    d = os.path.dirname(__file__)
    tf_df = pd.read_csv('{}/../refs/biomart_tf_gids.tsv'.format(d), sep='\t')
    tf_gids = tf_df['Gene stable ID'].unique().tolist()
    df['gid_stable'] = df['gid'].str.split('.', expand=True)[0]
    df.loc[df.gid_stable.isin(tf_gids), 'tf'] = True
    df.drop('gid_stable', axis=1, inplace=True)

    # and save
    df.to_csv(o, sep='\t', index=False)

def get_metadata_from_ab(df):
    """
    Get metadata dataframe from dagettaset names in
    TALON abundance file

    Parameters:
        df (pandas DataFrame): df of TALON abundance file
    Returns:
        df (pandas DataFrame): metadata for datasets

    """

    meta = pd.DataFrame()

    dataset_cols = get_dataset_cols(df)
    df = df[dataset_cols]
    df = df.transpose()
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)
    df = df['dataset'].to_frame()

    print(len(df.index))

    # label and add metadata for entries from
    # the mouse timecourse
    timepts = ['4d', '10d', '14d', '25d',
               '36d', '2mo', '18-20mo']
    c = 'b6cast'
    df['temp'] = df.dataset.str.split('_', expand=True)[1]
    df[c] = df['temp'].isin(timepts)
    df.drop('temp', axis=1, inplace=True)
    temp = df.loc[df[c]==True].copy(deep=True)
    temp[['tissue', 'age', 'sex', 'rep']] = temp['dataset'].str.split('_', expand=True)
    meta = pd.concat([meta, temp])

    # label and add metadata for entries from 5x vs wt brain
    df['temp'] = df.dataset.str.split('_', expand=True)[0]
    c = '5x_v_wt'
    df[c] = (df.b6cast==False)&\
            (df['temp'].isin(['hippocampus', 'cortex']))
    df.drop('temp', axis=1, inplace=True)
    temp = df.loc[df[c] == True].copy(deep=True)
    temp[['tissue', 'disease_status', 'sex']] = temp['dataset'].str.split('_', expand=True)[[0,1,2]]
    meta = pd.concat([meta, temp])

    # forelimb
    df['temp'] = df.dataset.str.split('_', expand=True)[0]
    c = 'forelimb'
    df[c] = (df.b6cast==False)&\
            (df['temp']=='forelimb')
    df.drop('temp', axis=1, inplace=True)
    temp = df.loc[df[c] == True].copy(deep=True)
    temp[['tissue', 'age']] = temp['dataset'].str.split('_', expand=True)[[0,1]]
    meta = pd.concat([meta, temp])

    # c2c12
    df['temp'] = df.dataset.str.split('_', expand=True)[0]
    c = 'c2c12'
    df[c] = (df.b6cast==False)&\
            (df['temp']=='c2c12')
    df.drop('temp', axis=1, inplace=True)
    temp = df.loc[df[c] == True].copy(deep=True)
    temp['tissue'] = temp['dataset'].str.rsplit('_', n=2, expand=True)[0]
    meta = pd.concat([meta, temp])

    # adrenal
    df['temp'] = df.dataset.str.split('_', expand=True)[0]
    c = 'adrenal'
    df[c] = (df.b6cast==False)&\
            (df['temp']=='adrenal')
    df.drop('temp', axis=1, inplace=True)
    temp = df.loc[df[c] == True].copy(deep=True)
    temp['tissue'] = temp['dataset'].str.split('_', expand=True)[0]
    meta = pd.concat([meta, temp])

    print(len(meta.index))


    return meta

def get_domains_per_tid(df):
    pass

def reformat_hmmer(hmmfile, ref=False):
    infile = open(hmmfile, 'r')
    outfile = open('temp_hmmer', 'w')

    # write header to file first
    header = ['tid', 'accession', 'bias', 'bitscore', 'description',\
             'cluster_num', 'domain_exp_num', 'domain_reported_num',\
             'env_num', 'evalue', 'id', 'overlap_num', 'region_num']
    if ref:
        header[0] = 'pid'
    outfile.write('\t'.join(header)+'\n')

    for record in SearchIO.parse(infile, 'hmmer3-tab'):

        if not ref:
            tid = record.id.split(';')[1]
        else:
            tid = record.id

        for hit in record.hits:
            new_line = []
            new_line.append(tid)
            new_line.append(str(hit.accession))
            new_line.append(str(hit.bias))
            new_line.append(str(hit.bitscore))
            new_line.append(str(hit.description))
            new_line.append(str(hit.cluster_num))
            new_line.append(str(hit.domain_exp_num))
            new_line.append(str(hit.domain_reported_num))
            new_line.append(str(hit.env_num))
            new_line.append(str(hit.evalue))
            new_line.append(str(hit.id))
            new_line.append(str(hit.overlap_num))
            new_line.append(str(hit.region_num))
            outfile.write('\t'.join(new_line)+'\n')
    infile.close()
    outfile.close()

    df = pd.read_csv('temp_hmmer', sep='\t')
    os.remove('temp_hmmer')

    # merge with tid/gid info
    f = '/Users/fairliereese/mortazavi_lab/data/mousewg/refs/gencode_vM21_pid_to_tid.tsv'
    tids = pd.read_csv(f, sep='\t')
    if ref:
        df = df.merge(tids, on='pid')
    else:
        df = df.merge(tids, on='tid')
    return df

def get_dataset_cols(df,
                    sample=None):
    """
    Get the names of the dataset columns from a TALON abundance file

    Parameters:
        df (pandas DataFrame): TALON ab dataframe
        sample (str): Tissue name of samples to pass
    Returns:
        dataset_cols (list of str): Names of dataset columns
    """

    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                           'annot_transcript_id', 'annot_gene_name',
                           'annot_transcript_name', 'n_exons', 'length',
                           'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]

    # subset by sample if necessary
    if sample != 'all':
        if sample == 'mouse_match':
            pass
        elif sample:
            dataset_cols = [x for x in dataset_cols if sample in x]
    return dataset_cols

def get_lr_bulk_abundance(datasets, filt=False, opref=None):
    """
    Given a list of datasets, return the TALON abundance table
    restricted to that list of datasets.

    Parameters:
        datasets (list of str): Datasets to include {'cortex', 'hippocampus',
            'heart', 'gastroc', 'adrenal'}
        filt (bool): Whether to return filtered or unfiltered.
            Default: False
        opref (str): Output file prefix to save as
            Default: None
    Returns:
        df (pandas DataFrame): TALON ab table with only those datasets included
    """

    d = os.path.dirname(__file__)
    if filt:
        fname = '{}/../lr_bulk/talon/mouse_talon_abundance_filtered.tsv'.format(d)
    else:
        fname = '{}/../lr_bulk/talon/mouse_talon_abundance.tsv'.format(d)

    df = pd.read_csv(fname, sep='\t')

    non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                           'annot_transcript_id', 'annot_gene_name',
                           'annot_transcript_name', 'n_exons', 'length',
                           'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_columns ]
    cols = []
    for d in datasets:
        temp = [c for c in dataset_cols if d in c]
        cols += temp
    cols = list(set(cols))
    cols.sort()
    cols = non_dataset_columns+cols

    df = df[cols]

    if opref:
        if filt:
            fname = '{}_talon_abundance_filtered.tsv'.format(opref)
        else:
            fname = '{}_talon_abundance.tsv'.format(opref)
        df.to_csv(fname, sep='\t', index=False)
    return df

def assign_gisx_sector(df):
    if 'n_tss' in df.columns.tolist():
        tss = 'n_tss'
        tes = 'n_tes'
    else:
        tss = 'tss'
        tes = 'tes'
    spl = 'splicing_ratio'
    df['total'] = df[tss]+df[tes]+df[spl]
    df['tss_ratio'] = df[tss] / df.total
    df['tes_ratio'] = df[tes] / df.total
    df['spl_ratio'] = df[spl] / df.total

    df['sector'] = 'simple'

    df.loc[df.tss_ratio > 0.5, 'sector'] = 'tss'
    df.loc[df.tes_ratio > 0.5, 'sector'] = 'tes'
    df.loc[df.spl_ratio > 0.5, 'sector'] = 'splicing'

    # mixed genes
    df.loc[(df.sector=='simple')&(df.n_iso>1), 'sector'] = 'mixed'

    return df


def filter_cells(adata, min_umi,
                 max_umi,
                 max_mt,
                 min_genes,
                 depth, verbose=False):
    """
    """
    if verbose:
        n = len(adata.obs.loc[adata.obs.depth == depth].index)
        print('# cells for depth {}: {}'.format(depth, n))

    # either has to satisfy the cutoff or be at a different depth

    # > ### UMI filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.total_counts > min_umi)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
        n = len(adata.obs.loc[adata.obs.depth == depth].index)
        print('# cells for depth {} after removing cells with < {} UMI: {}'.format(depth, min_umi, n))

    # < ### UMI filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.total_counts < max_umi)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
       n = len(adata.obs.loc[adata.obs.depth == depth].index)
       print('# cells for depth {} after removing cells with > {} UMI: {}'.format(depth, max_umi, n))

    # > ### genes filter
    if min_genes:
        depth_inds = (adata.obs.depth != depth)
        inds = (adata.obs.n_genes_by_counts > min_genes)|(depth_inds)
        adata = adata[inds, :]
        if verbose:
           n = len(adata.obs.loc[adata.obs.depth == depth].index)
           print('# cells for depth {} after removing cells with < {} genes: {}'.format(depth, min_genes, n))

    # < % MT filter
    depth_inds = (adata.obs.depth != depth)
    inds = (adata.obs.pct_counts_mt < max_mt)|(depth_inds)
    adata = adata[inds, :]
    if verbose:
       n = len(adata.obs.loc[adata.obs.depth == depth].index)
       print('# cells for depth {} after removing cells with < {} % MT: {}'.format(depth, max_mt, n))

    return adata

def read_raw_data(meta, genes, mtx, depth, experiment, dataset):
    """
    Read raw data into an anndata and furnish metadata with given
    variables

    Parameters:
        meta (str): Path to metadata.csv file
        genes (str): Path to genes.csv file
        mtx (str): Path to expression.mtx file
        depth (str): Depth of sequencing for this experiment, 12k or 2k
        experiment (str): Experiment, normal, 8mo, or ont
        dataset (str): 'cortex', 'hippocampus', '8mo_cortex', '8mo_hippocampus',
            or 'ont_hippocampus'
    """

    obs = pd.read_csv(meta)
    var = pd.read_csv(genes)
    adata = sc.read_mtx(mtx)
    X = adata.X
    adata = anndata.AnnData(X=X, obs=obs, var=var)

    adata.obs['depth'] = depth
    adata.obs['experiment'] = experiment
    adata.obs.set_index('cell_barcode', inplace=True)
    adata.var.set_index('gene_id', inplace=True)

    # add additional metadata
    if '8mo' not in dataset:
        if dataset == 'ont_hippocampus':
            adata.obs['tissue'] = 'hippocampus'
        else:
            adata.obs['tissue'] = dataset
        adata.obs['age'] = adata.obs['sample'].str.split('_', expand=True)[1]
        adata.obs['sex'] = adata.obs['sample'].str.split('_', expand=True)[2]
        adata.obs['rep'] = adata.obs['sample'].str.split('_', expand=True)[3]
    else:
        adata.obs['tissue'] = dataset.split('_')[1]
        adata.obs['age'] = '8m'
        adata.obs['genotype'] = adata.obs['sample'].str.split('_', expand=True)[0]
        adata.obs['tissue_from_sample'] = adata.obs['sample'].str.split('_', expand=True)[1]
        sex_map = {'C_27': 'F',
                   'C_28': 'F',
                   'C_29': 'M',
                   'C_30': 'M',
                   'HC_27': 'F',
                   'HC_28': 'F',
                   'HC_29': 'M',
                   'HC_30': 'M',
                   'C_33': 'F',
                   'C_34': 'F',
                   'C_35': 'M',
                   'C_36': 'M',
                   'HC_33': 'F',
                   'HC_34': 'F',
                   'HC_35': 'M',
                   'HC_36': 'M'}
        rep_map = {'C_27': '1',
                   'C_28': '2',
                   'C_29': '1',
                   'C_30': '2',
                   'HC_27': '1',
                   'HC_28': '2',
                   'HC_29': '1',
                   'HC_30': '2',
                   'C_33': '1',
                   'C_34': '2',
                   'C_35': '1',
                   'C_36': '2',
                   'HC_33': '1',
                   'HC_34': '2',
                   'HC_35': '1',
                   'HC_36': '2'}
        adata.obs['tissue_rep'] = adata.obs['sample'].str.split('_', expand=True, n=1)[1]
        adata.obs['sex'] = adata.obs.tissue_rep.map(sex_map)
        adata.obs['rep'] = adata.obs.tissue_rep.map(rep_map)
        adata.obs['sample'] = adata.obs.tissue_from_sample+'_'+\
                              adata.obs.age+'_'+\
                              adata.obs.sex+'_'+\
                              adata.obs.rep

        # limit to valid samples
        adata = adata[adata.obs.loc[adata.obs.genotype == 'B'].index, :]
        if dataset == '8mo_cortex':
            adata = adata[adata.obs.loc[adata.obs.tissue_from_sample == 'C'].index, :]
        elif dataset == '8mo_hippocampus':
            adata = adata[adata.obs.loc[adata.obs.tissue_from_sample == 'HC'].index, :]

        # drop extra columns
        adata.obs.drop(['genotype', 'tissue_from_sample', 'tissue_rep'],
            axis=1,
            inplace=True)


    return adata

def format_cellbender_data(adata, dataset, depth):
    """
    Format cellbender data to use gene names and report depth

    Parameters:
        adata (anndata AnnData): the anndata`
        dataset (str): cortex, hippocampus, gastroc, adrenal, or heart
        depth (str): N input cells (ie 1k, 2k, 12k...)
    """
    if dataset == 'cortex':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/cortex/sr_splitseq/splitpipe/'
        if depth == '12k':
            genes = d+'cortex_12k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'cortex_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        elif depth == '2k':
            genes = d+'cortex_2k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'cortex_2k/all-well/DGE_unfiltered/cell_metadata.csv'

    elif dataset == 'hippocampus':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/splitpipe/'
        if depth == '12k':
            genes = d+'hippocampus_12k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'hippocampus_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        elif depth == '2k':
            genes = d+'hippocampus_2k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'hippocampus_2k/all-well/DGE_unfiltered/cell_metadata.csv'

    elif dataset == '8mo_cortex' or dataset == '8mo_hippocampus':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/8mo/sr_splitseq/splitpipe/'
        if depth == '12k':
            genes = d+'8mo_12k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'8mo_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        elif depth == '2k':
            genes = d+'8mo_2k/all-well/DGE_unfiltered/genes.csv'
            meta = d+'8mo_2k/all-well/DGE_unfiltered/cell_metadata.csv'

    elif dataset == 'ont_hippocampus':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/splitpipe/'
        genes = d+'hippocampus_ont_match_2k/all-well/DGE_unfiltered/genes.csv'
        meta = d+'hippocampus_ont_match_2k/all-well/DGE_unfiltered/cell_metadata.csv'

    # add real gene names
    adata.var['gene_id'] = adata.var.index
    gene_df = pd.read_csv(genes)
    gene_df = gene_df[['gene_id', 'gene_name']]

    # merge with anndata var
    temp = adata.var.merge(gene_df, how='left', on='gene_id')
    temp.index = temp.gene_name
    adata.var = temp
    adata.var.index.name = None
    adata.var.drop('gene_ids', axis=1, inplace=True)
    adata.var_names_make_unique()

    # add sample / cell metadata
    df = pd.read_csv(meta)
    df = df[['cell_barcode', 'sample', 'rnd1_well', 'rnd2_well', 'rnd3_well']]
    adata.obs['cell_barcode'] = adata.obs.index

    # merge with anndata obs
    temp = adata.obs.merge(df, how='left', on='cell_barcode')
    temp.index = temp.cell_barcode
    adata.obs = temp
    adata.obs.index.name = None

    # add additional metadata

    if '8mo' not in dataset:
        if dataset == 'ont_hippocampus':
            adata.obs['tissue'] = 'hippocampus'

            # fix sample swaps
            # HC_10_F_2 swap with HC_10_M_2 from mine and Liz's conversation
            print(adata.obs.loc[adata.obs['sample']=='HC_10_F_2'].head())
            print(adata.obs.loc[adata.obs['sample']=='HC_10_M_2'].head())

            adata.obs.loc[adata.obs['sample']=='HC_10_F_2', 'sample'] = 'HC_10_M_2_temp'
            adata.obs.loc[adata.obs['sample']=='HC_10_M_2', 'sample'] = 'HC_10_F_2_temp'
            adata.obs.loc[adata.obs['sample']=='HC_10_F_2_temp', 'sample'] = 'HC_10_F_2'
            adata.obs.loc[adata.obs['sample']=='HC_10_M_2_temp', 'sample'] = 'HC_10_M_2'

            print(adata.obs.loc[adata.obs['sample']=='HC_10_F_2'].head())
            print(adata.obs.loc[adata.obs['sample']=='HC_10_M_2'].head())

        else:
            adata.obs['tissue'] = dataset
        adata.obs['age'] = adata.obs['sample'].str.split('_', expand=True)[1]
        adata.obs['sex'] = adata.obs['sample'].str.split('_', expand=True)[2]
        adata.obs['rep'] = adata.obs['sample'].str.split('_', expand=True)[3]

    # 8mo data has different sample structure than everything else
    else:
        adata.obs['tissue'] = dataset.split('_')[1]
        adata.obs['age'] = '8m'
        adata.obs['genotype'] = adata.obs['sample'].str.split('_', expand=True)[0]
        adata.obs['tissue_from_sample'] = adata.obs['sample'].str.split('_', expand=True)[1]
        sex_map = {'C_27': 'F',
                   'C_28': 'F',
                   'C_29': 'M',
                   'C_30': 'M',
                   'HC_27': 'F',
                   'HC_28': 'F',
                   'HC_29': 'M',
                   'HC_30': 'M',
                   'C_33': 'F',
                   'C_34': 'F',
                   'C_35': 'M',
                   'C_36': 'M',
                   'HC_33': 'F',
                   'HC_34': 'F',
                   'HC_35': 'M',
                   'HC_36': 'M'}
        rep_map = {'C_27': '1',
                   'C_28': '2',
                   'C_29': '1',
                   'C_30': '2',
                   'HC_27': '1',
                   'HC_28': '2',
                   'HC_29': '1',
                   'HC_30': '2',
                   'C_33': '1',
                   'C_34': '2',
                   'C_35': '1',
                   'C_36': '2',
                   'HC_33': '1',
                   'HC_34': '2',
                   'HC_35': '1',
                   'HC_36': '2'}
        adata.obs['tissue_rep'] = adata.obs['sample'].str.split('_', expand=True, n=1)[1]
        adata.obs['sex'] = adata.obs.tissue_rep.map(sex_map)
        adata.obs['rep'] = adata.obs.tissue_rep.map(rep_map)
        adata.obs['sample'] = adata.obs.tissue_from_sample+'_'+\
                              adata.obs.age+'_'+\
                              adata.obs.sex+'_'+\
                              adata.obs.rep

        # limit to valid samples
        adata = adata[adata.obs.loc[adata.obs.genotype == 'B'].index, :]
        if dataset == '8mo_cortex':
            adata = adata[adata.obs.loc[adata.obs.tissue_from_sample == 'C'].index, :]
        elif dataset == '8mo_hippocampus':
            adata = adata[adata.obs.loc[adata.obs.tissue_from_sample == 'HC'].index, :]

        # drop extra columns
        adata.obs.drop(['genotype', 'tissue_from_sample', 'tissue_rep'],
            axis=1,
            inplace=True)

    adata.obs['depth'] = depth

    # calc some metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = adata.var.gene_name.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    return adata

def load_cellbender_data(dataset):
    """
    Load / format all cellbender h5ad file metadata
    """
    if dataset == 'brain':
        hc = load_cellbender_data('hippocampus')
        ctx = load_cellbender_data('cortex')
        adata = hc.concatenate([ctx])

    elif dataset == 'cortex':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/cortex/sr_splitseq/cellbender/'
        experiment = 'normal'

        # shallow
        adata_shallow = sc.read(d+'12k/cortex_cellbender.h5ad')
        depth = '12k'
        adata_shallow = format_cellbender_data(adata_shallow, dataset, depth)
        adata_shallow.obs['experiment'] = experiment

        # deep
        adata_deep = sc.read(d+'2k/cortex_cellbender.h5ad')
        depth = '2k'
        adata_deep = format_cellbender_data(adata_deep, dataset, depth)
        adata_deep.obs['experiment'] = experiment

        # 8mo
        dataset_8mo = '8mo_cortex'
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/8mo/sr_splitseq/cellbender/'
        experiment = '8mo'

        # shallow
        adata_8mo_shallow = sc.read(d+'12k/8mo_cellbender.h5ad')
        depth = '12k'
        adata_8mo_shallow = format_cellbender_data(adata_8mo_shallow,
            dataset_8mo, depth)
        adata_8mo_shallow.obs['experiment'] = experiment

        # deep
        adata_8mo_deep = sc.read(d+'2k/8mo_cellbender.h5ad')
        depth = '2k'
        adata_8mo_deep = format_cellbender_data(adata_8mo_deep,
            dataset_8mo, depth)
        adata_8mo_deep.obs['experiment'] = experiment

        # concatenate
        adata = adata_shallow.concatenate([adata_deep,
                                           adata_8mo_deep,
                                           adata_8mo_shallow])

    elif dataset == 'hippocampus':
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/cellbender/'
        experiment = 'normal'

        # shallow
        adata_shallow = sc.read(d+'12k/hippocampus_cellbender.h5ad')
        depth = '12k'
        adata_shallow = format_cellbender_data(adata_shallow, dataset, depth)
        adata_shallow.obs['experiment'] = experiment

        # deep
        adata_deep = sc.read(d+'2k/hippocampus_cellbender.h5ad')
        depth = '2k'
        adata_deep = format_cellbender_data(adata_deep, dataset, depth)
        adata_deep.obs['experiment'] = experiment

        # 8mo
        dataset_8mo = '8mo_hippocampus'
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/8mo/sr_splitseq/cellbender/'
        experiment = '8mo'

        # shallow
        adata_8mo_shallow = sc.read(d+'12k/8mo_cellbender.h5ad')
        depth = '12k'
        adata_8mo_shallow = format_cellbender_data(adata_8mo_shallow,
            dataset_8mo, depth)
        adata_8mo_shallow.obs['experiment'] = experiment

        # deep
        adata_8mo_deep = sc.read(d+'2k/8mo_cellbender.h5ad')
        depth = '2k'
        adata_8mo_deep = format_cellbender_data(adata_8mo_deep,
            dataset_8mo, depth)
        adata_8mo_deep.obs['experiment'] = experiment

        # ONT matched 2k
        dataset_ont = 'ont_hippocampus'
        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/cellbender/'
        experiment = 'ont'

        adata_ont = sc.read(d+'ont_match_2k/hippocampus_cellbender.h5ad')
        depth = '2k'
        adata_ont = format_cellbender_data(adata_ont,
            dataset_ont, depth)
        adata_ont.obs['experiment'] = experiment

        # concatenate
        adata = adata_shallow.concatenate([adata_deep,
                                           adata_8mo_deep,
                                           adata_8mo_shallow,
                                           adata_ont])
        # adata = adata_shallow.concatenate([adata_deep,
        #                                    adata_8mo_deep,
        #                                    adata_8mo_shallow])

    return adata

def load_raw_data(dataset):
    """
    """
    # different files for each of the tissues
    # TODO - can probably automate this a lil more
    if dataset == 'brain':
        hc_adata = load_raw_data('hippocampus')
        ctx_adata = load_raw_data('cortex')
        adata = hc_adata.concatenate(ctx_adata)

    elif dataset == 'cortex':

        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/cortex/sr_splitseq/splitpipe/'

        # shallow
        meta = d+'cortex_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_12k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_12k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        experiment = 'normal'
        dset = 'cortex'
        adata_shallow = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # deep
        meta = d+'cortex_2k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_2k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_2k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        experiment = 'normal'
        dset = 'cortex'
        adata_deep = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # 8mo shallow
        meta = d+'cortex_12k_8mo/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_12k_8mo/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_12k_8mo/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        experiment = '8mo'
        dset = '8mo_cortex'
        adata_8mo_shallow = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # 8mo deep
        meta = d+'cortex_2k_8mo/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'cortex_2k_8mo/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'cortex_2k_8mo/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        experiment = '8mo'
        dset = '8mo_cortex'
        adata_8mo_deep = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        adata = adata_shallow.concatenate(adata_deep)
        adata = adata.concatenate(adata_8mo_shallow)
        adata = adata.concatenate(adata_8mo_deep)

        # # drop some trash
        # adata.obs.drop(['batch', 'Unnamed: 0'], axis=1, inplace=True)
        # adata.var.drop(['Unnamed: 0-0', 'Unnamed: 0-1'], axis=1, inplace=True)

    elif dataset == 'hippocampus':

        d = '/dfs6/pub/freese/mortazavi_lab/data/mousewg/hippocampus/sr_splitseq/splitpipe/'

        # shallow
        meta = d+'hippocampus_12k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_12k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_12k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        experiment = 'normal'
        dset = 'hippocampus'
        adata_shallow = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # deep
        meta = d+'hippocampus_2k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_2k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_2k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        experiment = 'normal'
        dset = 'hippocampus'
        adata_deep = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # 8mo shallow
        meta = d+'hippocampus_12k_8mo/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_12k_8mo/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_12k_8mo/all-well/DGE_unfiltered/DGE.mtx'
        depth = '12k'
        experiment = '8mo'
        dset = '8mo_hippocampus'
        adata_8mo_shallow = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # 8mo deep
        meta = d+'hippocampus_2k_8mo/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_2k_8mo/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_2k_8mo/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        experiment = '8mo'
        dset = '8mo_hippocampus'
        adata_8mo_deep = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        # ont match deep
        meta = d+'hippocampus_ont_match_2k/all-well/DGE_unfiltered/cell_metadata.csv'
        genes = d+'hippocampus_ont_match_2k/all-well/DGE_unfiltered/genes.csv'
        mtx = d+'hippocampus_ont_match_2k/all-well/DGE_unfiltered/DGE.mtx'
        depth = '2k'
        experiment = 'ont'
        dset = 'ont_hippocampus'
        adata_ont = read_raw_data(meta, genes, mtx, depth, experiment, dset)

        adata = adata_shallow.concatenate(adata_deep)
        adata = adata.concatenate(adata_8mo_shallow)
        adata = adata.concatenate(adata_8mo_deep)
        adata = adata.concatenate(adata_ont)

        # # drop some trash
        # adata.obs.drop(['batch', 'Unnamed: 0'], axis=1, inplace=True)
        # adata.var.drop(['Unnamed: 0-0', 'Unnamed: 0-1'], axis=1, inplace=True)

    # # add metadata
    # # add additional metadata
    # adata.obs['tissue'] = dataset
    # adata.obs['age'] = adata.obs['sample'].str.split('_', expand=True)[1]
    # adata.obs['sex'] = adata.obs['sample'].str.split('_', expand=True)[2]
    # adata.obs['rep'] = adata.obs['sample'].str.split('_', expand=True)[3]

    # calc some metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = adata.var.gene_name.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    return adata

def get_reads_per_bc(df):
    """
    Parameters:
        df (pandas DataFrame): read_annot df
    """
    df = df[['read_name', 'dataset']]
    temp = df.groupby('dataset').count().reset_index()
    temp.rename({'read_name': 'counts',
                 'dataset': 'barcode'}, axis=1, inplace=True)
    temp.sort_values(by='counts', ascending=False, inplace=True)
    return temp

def get_transcript_exp(df,
                       nov=['Known', 'NIC', 'NNC'],
                       filt_sc_novel=True,
                       filt_bcs=None,
                       filt_reads=None):
    """
    Parameters:
        df (pandas DataFrame): talon ab

    Parameters:
        df (pandas DataFrame): talon ab or read annot
        nov (list of str): List of novelty categories to include
        filt_sc_novel (bool): Filter out novel transcripts that are
            unique to the sc data (ie just include novel genes from
            the bulk data)
        filt_bcs (str): path to *_filt_bcs.txt file to filter on from
            read_annot file
    """
    # read annot
    if 'dataset' in df.columns:

        # remove novel genes
        df = df.loc[~df.annot_gene_id.str.contains('lr_splitseq')]

        if nov:
            df = df.loc[df.transcript_novelty.isin(nov)]

        if filt_sc_novel:
            df = df.loc[~df.annot_transcript_id.str.contains('lr_splitseq')]

        # filter out bcs
        if filt_bcs:
            bcs = pd.read_csv(filt_bcs, sep='\t', header=None)
            bcs.columns = ['bc']
            bcs = bcs.bc.tolist()
            df = df.loc[df.dataset.isin(bcs)]

        # filter out reads
        if filt_reads:
            df = df.loc[df.read_name.isin(filt_reads)]

        # limit to relevant columns
        cols = ['annot_transcript_id', 'annot_transcript_name', 'transcript_novelty',
                'annot_gene_name', 'annot_gene_id',
                'dataset', 'read_name']
        df = df[cols]

        # groupby dataset and count occurrences
        gb_cols = ['dataset', 'annot_transcript_id', 'annot_transcript_name', 'transcript_novelty',
                   'annot_gene_id', 'annot_gene_name']
        df = df.groupby(gb_cols).count()
        df.rename({'read_name': 'counts'}, axis=1, inplace=True)

        df.reset_index(inplace=True)
        df = df.pivot(index=['annot_transcript_id', 'annot_transcript_name',
                             'transcript_novelty', 'annot_gene_id', 'annot_gene_name'],
                 columns=['dataset'],
                 values=['counts'])
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''
        df.reset_index(inplace=True)
        df.fillna(0, inplace=True)
        df.reset_index(inplace=True)
        df.drop('index', inplace=True, axis=1)
        t_df = df

#     # talon abundance file
#     else:
#         # get gene level
#         non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
#                                'annot_transcript_id', 'annot_gene_name',
#                                'annot_transcript_name', 'n_exons', 'length',
#                                'gene_novelty', 'transcript_novelty', 'ISM_subtype']
#         dataset_cols = [ x for x in list(df.columns) \
#                             if x not in non_dataset_columns ]
#         id_col = ['annot_gene_id', 'annot_gene_name']
#         novelty_col = 'gene_novelty'

#         if filter_novel:
#             # only use known genes
#             df = df.loc[df.gene_novelty == 'Known']

    return t_df


def get_gene_exp(df,
                 filter_novel=True,
                 filt_bcs=None,
                 filt_reads=None):
    """
    Parameters:
        df (pandas DataFrame): talon ab or read annot
        filter_novel (bool): whether or not to filter out novel genes
        filt_bcs (str): path to *_filt_bcs.txt file to filter on from
            read_annot file
    """
    # read annot
    if 'dataset' in df.columns:

        if filter_novel:
            df = df.loc[df.gene_novelty == 'Known']

        if filt_bcs:
            bcs = pd.read_csv(filt_bcs, sep='\t', header=None)
            bcs.columns = ['bc']
            bcs = bcs.bc.tolist()
            df = df.loc[df.dataset.isin(bcs)]

        if filt_reads:
            df = df.loc[df.read_name.isin(filt_reads)]

        # limit to relevant columns
        cols = ['annot_gene_id', 'annot_gene_name', 'dataset', 'read_name']
        df = df[cols]

        # groupby gene name, id, and dataset and count occurrences
        gb_cols = ['annot_gene_id', 'annot_gene_name', 'dataset']
        df = df.groupby(gb_cols).count()
        df.rename({'read_name': 'counts'}, axis=1, inplace=True)

        df.reset_index(inplace=True)
        df = df.pivot(index=['annot_gene_id', 'annot_gene_name'],
                 columns=['dataset'],
                 values=['counts'])
        df.columns = df.columns.get_level_values(1)
        df.columns.name = ''
        df.reset_index(inplace=True)
        df.fillna(0, inplace=True)
        df.reset_index(inplace=True)
        df.drop('index', inplace=True, axis=1)
        gene_df = df

    # talon abundance file
    else:
        # get gene level
        non_dataset_columns = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                               'annot_transcript_id', 'annot_gene_name',
                               'annot_transcript_name', 'n_exons', 'length',
                               'gene_novelty', 'transcript_novelty', 'ISM_subtype']
        dataset_cols = [ x for x in list(df.columns) \
                            if x not in non_dataset_columns ]
        id_col = ['annot_gene_id', 'annot_gene_name']
        novelty_col = 'gene_novelty'

        if filter_novel:
            # only use known genes
            df = df.loc[df.gene_novelty == 'Known']

        # aggregate by gene
        gene_df = df[id_col+dataset_cols].groupby(id_col).sum()
        gene_df.reset_index(inplace=True)

    return gene_df

def get_bc1_matches():
    # from spclass.py - barcodes and their well/primer type identity
    d = os.path.dirname(__file__)
    bc_file = '{}/../refs/bc_8nt_v2.csv'.format(d)
    bc_df = pd.read_csv(bc_file, index_col=0, names=['bc'])
    bc_df['well'] = [i for i in range(0, 48)]+[i for i in range(0, 48)]
    bc_df['primer_type'] = ['dt' for i in range(0, 48)]+['randhex' for i in range(0, 48)]

    # pivot on well to get df that matches bcs with one another from the same well
    bc_df = bc_df.pivot(index='well', columns='primer_type', values='bc')
    bc_df = bc_df.rename_axis(None, axis=1).reset_index()
    bc_df.rename({'dt': 'bc1_dt', 'randhex': 'bc1_randhex'}, axis=1, inplace=True)

    return bc_df

def add_bc_info_illumina(sr_df=None,
                         dataset=None,
                         which='raw'):
    """
    Get illumina metadata with barcode information added to it

    Parameters:
        sr_df (pandas DataFrame): Dataframe with column 'bc' which contains 16ntbc_well#
        dataset (str): Choose from adrenal, cortex, hippocampus
        which (str): Choose from raw or processed
    """
    # we want info about primer type as well
    bc_df = get_bc1_matches()

    # if given a df, just add 24nt bc to df
    if isinstance(sr_df, pd.DataFrame):

        sr_df['rnd1_well'] = sr_df.bc.str.split('_', expand=True)[1].astype(int)

        # get 24nt barcode
        sr_df = sr_df.merge(bc_df, how='left', \
                     left_on='rnd1_well', right_on='well')

        sr_df.rename({'bc':'cell_barcode'}, axis=1, inplace=True)

        sr_df['bc3'] = sr_df.cell_barcode.str.slice(start=0, stop=8)
        sr_df['bc2'] = sr_df.cell_barcode.str.slice(start=8, stop=16)
        sr_df['bc1'] = sr_df.bc1_dt
        sr_df['bc'] = sr_df.bc3+sr_df.bc2+sr_df.bc1

    elif which == 'raw':
        # read in illumina bcs
        d = os.path.dirname(__file__)
        fname = '{}/../{}/sr_splitseq/scanpy/illumina_raw_metadata.csv'.format(d, dataset)
        sr_df = pd.read_csv(fname)

        # remove these random columns???
        drop_cols = ['rnd3_rnd2','Index','Barcode','cell_barcode_24nt',
                     'star_barcode','well','bc1_dt','bc1_randhex',
                     'cell_barcode_24nt_dt','cell_barcode_24nt_rh',
                     'star_barcode_dt','star_barcode_rh']
        drop_cols = list(set(sr_df.columns.tolist())&set(drop_cols))
        sr_df.drop(drop_cols, axis=1, inplace=True)

        # get 24nt barcode
        sr_df = sr_df.merge(bc_df, how='left', \
                     left_on='rnd1_well', right_on='well')
        sr_df['bc3'] = sr_df.cell_barcode.str.slice(start=0, stop=8)
        sr_df['bc2'] = sr_df.cell_barcode.str.slice(start=8, stop=16)
        sr_df['bc1'] = sr_df.bc1_dt
        sr_df['bc'] = sr_df.bc3+sr_df.bc2+sr_df.bc1

        # get library
        m = {'ont_match': 'b', 'normal': 'a'}
        sr_df['library'] = sr_df.experiment.map(m)

    elif which == 'processed':
        print('u didnt implement this yet dummy')
        sr_df = None

    return sr_df

def get_sample_metadata(dataset):
    """
    Get bc1_dt, sample id, and sample metadata from the raw
    illumina metadata, which tells us which wells each sample
    was loaded into

    Parameters:
        dataset (str): Choose from adrenal, hippocampus, cortex
    """

    # use raw data for the sample-level metadata
    samp_df = add_bc_info_illumina(dataset=dataset, which='raw')

    # only need sample ID, bc1_dt, and library
    samp_df = samp_df[['sample', 'bc1_dt', 'library']].drop_duplicates()

    # fix adrenal sample names
    if samp_df['sample'].values[0].count('_') < 3:
        samp_df['temp1'] = samp_df['sample'].str[:-1]
        samp_df['temp2'] = samp_df['sample'].str[-1]
        samp_df['sample'] = samp_df.temp1+'_'+samp_df.temp2
        samp_df.drop(['temp1', 'temp2'], axis=1, inplace=True)

    # grab metadata from references
    d = os.path.dirname(__file__)
#     d = '/Users/fairliereese/mortazavi_lab/data/mousewg/scripts/'

    fname = '{}/../refs/age_metadata.tsv'.format(d)
    age = pd.read_csv(fname, sep='\t')

    fname = '{}/../refs/sex_metadata.tsv'.format(d)
    sex = pd.read_csv(fname, sep='\t')

    fname = '{}/../refs/tissue_metadata.tsv'.format(d)
    tissue = pd.read_csv(fname, sep='\t')

    samp_df[['tissue', 'age', 'sex', 'rep']] = samp_df['sample'].str.split('_', expand=True)

    samp_df = samp_df.merge(age, how='left', left_on='age', right_on='short')
    samp_df = samp_df.merge(tissue, how='left', left_on='tissue', right_on='short')

    samp_df = samp_df[['sample', 'tissue_desc', 'age_desc', 'sex', 'rep', 'bc1_dt', 'library']]
    samp_df.rename({'tissue_desc': 'tissue', 'age_desc': 'age'},
                  axis=1, inplace=True)

    return samp_df

def get_illumina_metadata(dataset, include_raw=True):
    """
    Get a table with metadata pertaining to each cell from the
    corresponding illumina

    Parameters:
        dataset (str): Choose from 'adrenal', 'hippocampus', 'cortex'

    Returns:
        sr_df (pandas DataFrame): DataFrame with 24nt barcode, umi and gene count,
            and cluster / celltype labels if available
    """

    # get raw umi and gene count from the unfiltered metadata
    if include_raw:
        raw_df = add_bc_info_illumina(dataset=dataset, which='raw')

        raw_df = raw_df[['bc', 'umi_count', 'gene_count']]
        raw_df.rename({'umi_count': 'sr_umi_count',
                       'gene_count': 'sr_gene_count'}, axis=1, inplace=True)

    # get annotation information from the filtered metadata
    d = os.path.dirname(__file__)
    fname = '{}/../sr_splitseq/all_tissues_C2C12_Parse_deep_TFs_mirhgs_chromreg_metadata.tsv'.format(d)
    proc_df = pd.read_csv(fname, sep='\t')

    # drop dupes for some reason
    proc_df = proc_df.drop_duplicates()

    # only parse
    proc_df = proc_df.loc[proc_df.technology == 'Parse']

    # only deep
    proc_df = proc_df.loc[proc_df.depth1 == 'deep']

    # get relevant tissue
    tissue_map = {'adrenal': 'Adrenal',
              'hippocampus': 'Hippocampus',
              'cortex': 'Cortex'}
    proc_df = proc_df.loc[proc_df.tissue == tissue_map[dataset]]

    # add a vs b library name
    m = {'normal': 'a', 'ont_match': 'b'}
    proc_df['library'] = proc_df.experiment.map(m)

    # rename bc column
    proc_df.rename({'cell_barcode_24nt': 'bc'}, axis=1, inplace=True)

    drop_cols = ['experiment', 'tissue', 'sex', 'timepoint',
                 'sample', 'doublets', 'rep']
    proc_df.drop(drop_cols, axis=1, inplace=True)

    # rename stuff
    proc_df.rename({'celltypes': 'sr_celltype',
                    'gen_celltype': 'sr_gen_celltype',
                    'subtypes': 'sr_sub_celltype'},
                   axis=1, inplace=True)

    if include_raw:
        sr_df = raw_df.merge(proc_df, how='left', on='bc')
    else:
        sr_df = proc_df.copy(deep=True)

    return sr_df

# def get_illumina_metadata(dataset, include_raw=True):
#     """
#     Get a table with metadata pertaining to each cell from the
#     corresponding illumina

#     Parameters:
#         dataset (str): Choose from 'adrenal', 'hippocampus', 'cortex'

#     Returns:
#         sr_df (pandas DataFrame): DataFrame with 24nt barcode, umi and gene count,
#             and cluster / celltype labels if available
#     """

#     # get raw umi and gene count from the unfiltered metadata
#     raw_df = add_bc_info_illumina(dataset=dataset, which='raw')

#     raw_df = raw_df[['bc', 'umi_count', 'gene_count']]
#     raw_df.rename({'umi_count': 'sr_umi_count',
#                    'gene_count': 'sr_gene_count'}, axis=1, inplace=True)

#     # get annotation information from the filtered metadata
#     d = os.path.dirname(__file__)
#     fname = '{}/../{}/sr_splitseq/scanpy/illumina_processed_metadata.tsv'.format(d, dataset)
#     proc_df = pd.read_csv(fname, sep='\t')

#     # each of the datasets has slightly different metadata
#     if dataset == 'adrenal':
#         proc_df.rename({'seurat_clusters': 'sr_clusters',
#                    'celltypes': 'sr_celltype',
#                    'gen_celltypes': 'sr_gen_celltype',
#                    'barcode': 'bc'}, axis=1, inplace=True)
# #         proc_df.drop('Sample', axis=1, inplace=True)
#     elif dataset == 'hippocampus':
#         proc_df = proc_df.loc[proc_df.tissue == 'Hippocampus']
#         proc_df['bc'] = proc_df['cellID'].str.split('.', n=1, expand=True)[0]
#         temp = proc_df.copy(deep=True)
#         temp.sort_values(by='bc')
#         proc_df.rename({'seurat_clusters': 'sr_clusters',
#                         'celltypes': 'sr_celltype',
#                         'gen_celltypes': 'sr_gen_celltype'},
#                          axis=1, inplace=True)
#         m = {'normal': 'a', 'ont_match': 'b'}
#         proc_df['library'] = proc_df.experiment.map(m)
#         drop_cols = ['experiment', 'tissue', 'sex', 'timepoint']
#         proc_df.drop(drop_cols, axis=1, inplace=True)
#         proc_df = add_bc_info_illumina(sr_df=proc_df)
#     elif dataset == 'cortex':
#         print('u havent done this yet')

#     print('raw counts + genes / cell will be incorrect until you replace metadata files')

#     if include_raw:
#         sr_df = raw_df.merge(proc_df, how='left', on='bc')
#     else:
#         sr_df = proc_df.copy(deep=True)

#     sr_df.sr_clusters = sr_df.sr_clusters.astype('category')
#     return sr_df

def add_obs_counts_col(adata, col):
    """
    Add a counts column and a counts string column for each category in
    an AnnData's obs
    """
    counts_col = col+'_counts'
    counts_str_col = col+'_counts_str'
    if counts_str_col in adata.obs.columns.tolist():
        print('Already added counts for {}'.format(col))
        return adata, None

    temp = adata.obs[['bc', col]].groupby(col).count().reset_index()
    temp.rename({'bc': counts_col}, axis=1, inplace=True)

    adata.obs['bc_index_back'] = adata.obs.index
    adata.obs = adata.obs.merge(temp, on=col, how='left')
    adata.obs.rename({'bc_index_back': 'bc_index'}, axis=1, inplace=True)
    adata.obs.set_index('bc_index', inplace=True)
    adata.obs[counts_str_col] = adata.obs[col].astype(str)+' ('+adata.obs[counts_col].astype(str)+')'
    adata.obs.loc[adata.obs[counts_str_col] == 'nan (nan)', counts_str_col] = np.nan
    return adata, temp

def make_adata(df,
               samp_df,
               sr_df,
               verbose=False,
               how='gene'):

    """
    Make an AnnData from long read data.

    Parameters:
        df (pandas DataFrame): DataFrame from get_gene_exp function
        samp_df (pandas DataFrame): DataFrame from get_sample_metadata function
        sr_df (pandas DataFrame): DataFrame from get_illumina_metadata function
        verbose (bool): Whether or not to display output messages
        how (str): Choose from 'gene' or 'transcript'
    """

    if how == 'gene':
        # print(df.head())
        var = df[['annot_gene_id', 'annot_gene_name']]
        df.drop(['annot_gene_name'], axis=1, inplace=True)
        df.set_index('annot_gene_id', inplace=True)
    elif how == 'transcript':
        var = df[['annot_transcript_id', 'annot_transcript_name', \
                'annot_gene_id', 'annot_gene_name', 'transcript_novelty']]
        df.drop(['annot_transcript_name', 'annot_gene_id', \
                 'annot_gene_name', 'transcript_novelty'], axis=1, inplace=True)
        df.set_index('annot_transcript_id', inplace=True)

    df = df.transpose()
    df.index.name = 'bc'
    X = scipy.sparse.csr_matrix(df.values)
    df.reset_index(inplace=True)
    obs = df.bc.to_frame()
    obs = df.bc.to_frame()
    obs['bc3_long'] = obs['bc'].str.slice(0,8)
    obs['bc2_long'] = obs['bc'].str.slice(8,16)
    temp = obs['bc'].values[0]
    if '-' in temp:
        obs['bc1_long'] = obs['bc'].str.slice(16,24)
        obs['bc_index'] = obs.bc
        obs['bc_back'] = obs.bc.str.split('-', expand=True)[0]
        obs['experiment'] = obs.bc.str.split('-', expand=True)[1]
        obs['library'] = obs.experiment.str.split('_', expand=True)[1]

        # fix library names for those with issues
        if '2ka' in obs.library.unique().tolist():
            obs.loc[obs.library == '2ka', 'library'] = 'a'

        obs['bc'] = obs.bc_back
        obs.drop('bc_back', axis=1, inplace=True)

    else:
        obs['bc1_long'] = obs['bc'].str.slice(16,-1)
        obs['bc_index'] = obs['bc']

    if verbose:
        print('Found {} unique bc3s'.format(len(obs.bc3_long.unique())))
        print('Found {} unique bc2s'.format(len(obs.bc2_long.unique())))
        print('Found {} unique bc1s'.format(len(obs.bc1_long.unique())))

    # merge with information from illumina runs
    # pdb.set_trace()

    if sr_df is not None:
        # print(len(obs.index))
        obs = obs.merge(sr_df, how='left', on=['bc', 'library'])
        # print(len(obs.index))

    # merge with sample information
    if samp_df is not None:
        # print()
        # print(len(obs.index))
        obs = obs.merge(samp_df, how='left', left_on=['bc1_long', 'library'], right_on=['bc1_dt', 'library'])
        # print(len(obs.index))

    # construct the actual anndata
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obs.set_index('bc_index', inplace=True)

    if how == 'gene':
        adata.var.set_index('annot_gene_id', inplace=True)
    elif how == 'transcript':
        adata.var.set_index('annot_transcript_id', inplace=True)
    adata.layers['raw'] = adata.X.copy()

    # annotate the group of mitochondrial genes as 'mt'
    adata.var['mt'] = adata.var.annot_gene_name.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # format gene names
    if how == 'gene':
        adata.var.reset_index(inplace=True)
        adata.var.annot_gene_name = adata.var.annot_gene_name.astype(str)
        adata.var.index = adata.var.annot_gene_name
        adata.var_names_make_unique()
        adata.var.drop('annot_gene_name', axis=1, inplace=True)
    elif how == 'transcript':
        adata.var.reset_index(inplace=True)
        adata.var.annot_transcript_name = adata.var.annot_transcript_name.astype(str)
        adata.var.index = adata.var.annot_transcript_name
        adata.var_names_make_unique()
        adata.var.drop('annot_transcript_name', axis=1, inplace=True)

    return adata


def write_seurat_tables(adata, opref):
    obs = adata.obs
    x = adata.X

    x_coords = [n[0] for n in adata.obsm['X_umap']]
    y_coords = [n[1] for n in adata.obsm['X_umap']]

    obs['umap_x'] = x_coords
    obs['umap_y'] = y_coords

    for i in range(adata.obsm['X_pca'].shape[1]):
        index = 'pca_' + str((i+1))
        coords = [n[i] for n in adata.obsm['X_pca']]
        obs[index] = x_coords

    row_labels = obs.index.tolist()
    col_labels = adata.var.index.tolist()

    x = pd.DataFrame(data=x, index=row_labels, columns=col_labels)

    obs.to_csv('{}_obs.tsv'.format(opref), sep='\t')
    x.to_csv('{}_x.tsv'.format(opref), sep='\t')

def write_swan_input_tables(adata, opref):

    # dump to formats that swan can handle
    meta = adata.obs.copy(deep=True)
    meta.reset_index(inplace=True)
    meta.rename({'bc':'dataset'}, axis=1, inplace=True)
    fname = '{}_swan_metadata.tsv'.format(opref)
    meta.to_csv(fname, sep='\t', index=False)

    index = adata.obs.index.tolist()
    cols = adata.var.index.tolist()
    data = adata.layers['raw']

    print(len(cols))
    print(len(index))
    df = pd.DataFrame(data=data, index=index, columns=cols)
    df.reset_index(inplace=True)
    df.rename({'index': 'dataset'}, axis=1, inplace=True)
    df.set_index('dataset', inplace=True)
    df = df.transpose()
    df.index.name = 'transcript_id'
    df.reset_index(inplace=True)
    fname = '{}_swan_abundance.tsv'.format(opref)
    df.to_csv(fname, sep='\t', index=False)

def make_sc_sg(adata, dataset='adrenal', save=None):
    """
    Make a SwanGraph from LR-Split-seq data.

    Parameters:
        adata (anndata AnnData): AnnData object to pull abundance/metadata from
        dataset (str): Choose from 'adrenal', 'hippocampus'
        save (bool): Whether to save SwanGraph
    """

    annot = '../../../refs/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf'

    d = os.path.dirname(__file__)
#     d ='/Users/fairliereese/mortazavi_lab/data/mousewg/scripts/'
    fname = '{}/../{}/lr_splitseq/scanpy/transcript_processed.h5ad'.format(d, dataset)

    # use the transcriptome from bulk
    gtf = '{}/../lr_bulk/talon/mouse_all_known_nic_nnc_talon.gtf'.format(d)

    # counts matrix
    bcs = adata.obs.index.tolist()
    X = adata.layers['raw'].toarray()
    tids = adata.var.annot_transcript_id.tolist()
    df = pd.DataFrame(data=X, columns=tids, index=bcs)
    df = df.transpose()
    df.index.name = 'transcript_id'
    ab = '{}/../{}/lr_splitseq/scanpy/temp_ab.tsv'.format(d, dataset)
    df.to_csv(ab, sep='\t')

    # metadata
    obs = adata.obs.copy(deep=True)
    cols = ['tissue', 'age', 'sex', 'sr_celltype',
            'sr_gen_celltype', 'sr_sub_celltype',
            'leiden', 'lr_celltype', 'lr_sub_celltype', 'lr_gen_celltype', 'umap_x', 'umap_y']
    cols = list(set(obs.columns.tolist())&set(cols))
    obs = obs[cols]
    obs.reset_index(inplace=True)
    obs.rename({'bc_index': 'dataset'}, axis=1, inplace=True)
    meta = '{}/../{}/lr_splitseq/scanpy/temp_meta.tsv'.format(d, dataset)
    obs.to_csv(meta, sep='\t')

    # make the SwanGraph
    sg = swan.SwanGraph(sc=True)
    sg.add_annotation(annot)
    sg.add_transcriptome(gtf, include_isms=False)

    sg.save_graph('{}/../{}/lr_splitseq/swan/swan'.format(d, dataset))
    sg.add_abundance(ab)
    sg.add_metadata(meta)
    os.remove(ab)
    os.remove(meta)

    # save swangraph
    swan_dir = '{}/../{}/lr_splitseq/swan/'.format(d, dataset)
    try:
        os.mkdir(swan_dir)
    except:
        pass
    sg_pref = swan_dir+'swan'
    sg.save_graph(sg_pref)

    return sg

def write_lr_sc_gene_metadata(adata):
    """
    Write a table from the gene-level LR-Split-seq anndata
    that includes the metadata that we want to port over
    to the transcript-level LR-Split-seq anndata

    Parameters:
        adata (anndata AnnData): Gene-level AnnData obj.
            from LR-Split-seq
    """

    cols = ['leiden', 'lr_celltype', 'lr_sub_celltype']
    cols = list(set(cols)&set(adata.obs.columns.tolist()))
    obs = adata.obs[cols]

    x_coords = [n[0] for n in adata.obsm['X_umap']]
    y_coords = [n[1] for n in adata.obsm['X_umap']]
    obs['umap_x'] = x_coords
    obs['umap_y'] = y_coords

    # write umap coords, clustering results, cell bc,
    # and celltype
    obs.to_csv('gene_metadata.tsv', sep='\t')

def write_bc_leiden(adata, opref):
    obs = adata.obs['leiden']
    obs.to_csv('{}_leiden.tsv'.format(opref), sep='\t')

def rm_sirv_ercc(df):
    """From TALON ab file"""
    df = df.loc[~df.annot_gene_id.str.contains('SIRV')]
    df.loc[~df.annot_gene_id.str.contains('ERCC-')]
    return df

def get_counts_table(df,
                    how='gene',
                    nov=None,
                    min_tpm=None,
                    gene_subset=None,
                    save=False):
    """
    Parameters:
        df (pandas DataFrame): TALON abundance table
        sample (str): Choose from 'cell_line', 'tissue', or None
        how (str): Choose from 'gene' or 'iso'
        nov (list of str): List of accepted novelty types (w/ how='iso')
        min_tpm (float): Keep only genes / isos that have at least one
            TPM >= the value across the libraries
        gene_subset (str): Choose from 'polya' or None
        save (bool): Whether or not to save the output matrix

    Returns:
        df (pandas DataFrame): Counts for gene or isoforms in the requested
            samples above the input detection threshold.
        ids (list of str): List of str indexing the table
    """
    print('Calculating {} counts'.format(how))

    dataset_cols = get_dataset_cols(df)
    df = rm_sirv_ercc(df)

#     # merge with information about the gene
#     gene_df, _, _ = get_gtf_info(how='gene')
#     gene_df = gene_df[['gid', 'biotype_category', 'tf']]
#     df = df.merge(gene_df, how='left', left_on='annot_gene_id', right_on='gid')

    # get indices that we'll need to subset on
    if how == 'gene':
        id_col = 'annot_gene_id'
        nov_col = 'gene_novelty'
        nov = ['Known']
    elif how == 'iso':
        id_col = 'annot_transcript_id'
        nov_col = 'transcript_novelty'

    # filter on novelty
    if nov:
        print('Subsetting for novelty categories {}'.format(nov))
        nov_inds = df.loc[df[nov_col].isin(nov), id_col].tolist()
    else:
        nov_inds = df[id_col].tolist()

    # filter on gene subset
    if gene_subset:
        print('Subsetting for {} genes'.format(gene_subset))
        if gene_subset == 'polya':
            polya_cats = ['protein_coding', 'pseudogene', 'lncRNA']
            gene_inds = df.loc[df.biotype_category.isin(polya_cats), id_col].tolist()
        elif gene_subset == 'tf':
            gene_inds = df.loc[df.tf == True, id_col].tolist()
    else:
        gene_inds = df[id_col].tolist()

    # get intersection of both
    subset_inds = list(set(nov_inds)&set(gene_inds))

    # sum up counts across the same gene
    if how == 'gene':
        df = df[dataset_cols+[id_col]]
        df = df.groupby(id_col).sum().reset_index()

    # set index so that all values in df reflect
    # counts per transcript or gene
    df.set_index(id_col, inplace=True)

#     # compute TPM
#     tpm_cols = []
#     for d in dataset_cols:
#         tpm_col = '{}_tpm'.format(d)
#         total_col = '{}_total'.format(d)
#         df[total_col] = df[d].sum()
#         df[tpm_col] = (df[d]*1000000)/df[total_col]
#         tpm_cols.append(tpm_col)
#     df = df[tpm_cols]

#     # reformat column names
#     df.columns = [c.rsplit('_', maxsplit=1)[0] for c in df.columns]

#     # enforce tpm threshold
#     if min_tpm:
#         print('Enforcing minimum TPM')
#         print('Total # {}s detected: {}'.format(how, len(df.index)))
#         df = df.loc[(df >= min_tpm).any(axis=1)]
#         print('# {}s >= {} tpm: {}'.format(how, min_tpm, len(df.index)))

    # subset if necessary
    if gene_subset or nov:
        print('Applying gene type and novelty subset')
        df = df.loc[df.index.isin(subset_inds)]

    print('Number of {}s reported: {}'.format(how, len(df.index)))

    if save:
        fname = '{}_{}_tpm.tsv'.format(sample, how)
        df.to_csv(fname, sep='\t')

    ids = df.index.tolist()

    return df, ids

def get_det_table(df,
                  how='gene',
                  min_tpm=1,
                  gene_subset='polya',
                  sample='all',
                  groupby='library',
                  nov='Known'):
    """
    Get a dataframe of True / False whether or not a gene / isoform
    was detected in a specific library or sample

    Parameters:
        df (pandas DataFrame): TALON abundance
        how (str): Either "gene" or "iso"
        min_tpm (float): Minimum TPM to call a gene / iso as detected
        gene_subset (str): Subset of genes to use, 'polya' or None
        groupby (str): Either "sample", 'library', or 'all'
            used to groupby datasets displayed
        nov (list of str): Only used with how='iso', novelty categories of
            isoforms to consider

    Returns:
        df (pandas DataFrame): DataFrame with True / False entries
            for each isoform / gene per library / sample
    """

    # calc TPM per library on desired samples
    df, tids = get_tpm_table(df,
                   how=how,
                   nov=nov,
                   min_tpm=min_tpm,
                   sample=sample,
                   gene_subset=gene_subset)

    df = df.transpose()
    df.index.name = 'dataset'
    df.reset_index(inplace=True)

    # set up df to groupby sample or library
    if groupby == 'sample':

        # add biosample name (ie without rep / sex information)
        df['biosample'] = df.dataset.str.rsplit('_', n=2, expand=True)[0]
        df.drop(['dataset'], axis=1, inplace=True)

        print('Found {} total samples'.format(len(df.biosample.unique().tolist())))
        df = df.groupby('biosample').max()

    elif groupby == 'library':
        df.rename({'dataset': 'library'}, axis=1, inplace=True)
        print('Found {} total libraries'.format(len(df.library.unique().tolist())))
        df = df.groupby('library').max()

    elif groupby == 'all':
        df['dataset'] = 'all'
        df = df.groupby('dataset').max()

    df = (df >= min_tpm)
    return df

def compute_triplets(t_df,
                     df,
                     min_tpm=1,
                     sample='all',
                     groupby='sample',
                     bool_col=None):
    """
    Compute the triplets on the sample or library level

    Parameters:
        t_df (pandas DataFrame): t_df output from get_ic_tss_tes
        df (pandas DataFrame): Filtered TALON abundance or 90% set
            dataframe
        min_tpm (int): Min TPM to be considered detected
        sample (str): Choose 'cell_line', 'tissue', 'mouse_match'
        groupby (str): Choose 'library', 'sample', or 'all'

    Returns:
        counts (pandas DataFrame): DF w/ n tss, ic, tes, and
            splicing ratio for each gene in each sample / lib
    """

    # get table of which isoforms are detected in
    # which samples / libraries, or which isoforms
    # are part of the 90% set / sample
    if 'transcript_novelty' in df.columns:
        df = get_det_table(df,
                           how='iso',
                           min_tpm=min_tpm,
                           sample=sample,
                           groupby=groupby,
                           nov=['Known', 'NIC', 'NNC'])

    # otherwise, expect 90% set file format and coerce into
    # boolean 90% set detection table format
    elif 'pi' in df.columns:
        df = df[['tid', 'biosample']]
        df['in_90_set'] = True
        df = df.pivot(index='biosample', columns='tid', values='in_90_set').fillna(value=False)
        df.columns.name = ''

    # loop through samples / libraries compute triplet
    # for detected transcripts in each sample
    counts = pd.DataFrame()
    for ind, entry in df.iterrows():
        entry = entry.to_frame()
        tids = entry.loc[entry[ind] == True].index.tolist()
        temp = count_tss_ic_tes(t_df, subset=tids)
        temp['source'] = ind
        if bool_col:
            temp[bool_col] = True
        counts = pd.concat([counts, temp])

    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')

    return counts

def get_gtf_info(how='gene',
                 subset=None,
                 add_stable_gid=False,
                 ver='v40_cerberus'):
    """
    Gets the info from the annotation about genes / transcripts

    Parameters:
        how (str): 'gene', 'iso', 'ic', 'tss', 'tes'
        subset (str): 'polya', 'tf', 'protein_coding' or None
        add_stable_gid (bool): Add stable gid (code from cerberus)
        ver (str): {'v29', 'v40_cerberus'}

    Returns:
        df (pandas DataFrame): DataFrame with info for gene / transcript
        biotype_counts (pandas DataFrame): DataFrame with the counts
            per biotype reported in gencode
        biotype_cat_counts (pandas DataFrame): DataFrame with the counts
            per meta biotype reported in gencode
    """
    iso_hows = ['iso', 'tss', 'tes', 'ic']
    d = os.path.dirname(__file__)
    if how == 'gene' and ver == 'v29':
        fname = '{}/../refs/gencode_v29_gene_metadata.tsv'.format(d)
    elif how in iso_hows and ver == 'v29':
        fname = '{}/../refs/gencode_v29_transcript_metadata.tsv'.format(d)
    elif how == 'gene' and ver == 'v40_cerberus':
        fname = '{}/../refs/cerberus/v40_gene_metadata.tsv'.format(d)
    elif how in iso_hows and ver == 'v40_cerberus':
        fname = '{}/../refs/cerberus/v40_transcript_metadata.tsv'.format(d)
    elif how == 'gene' and ver == 'vM21_cerberus':
        fname = '{}/../refs/gencode_vM21_gene_metadata.tsv'.format(d)
    elif how in iso_hows and ver == 'vM21_cerberus':
        fname = '{}/../refs/gencode_vM21_transcript_metadata.tsv'.format(d)
    elif how == 'gene' and ver == 'vM25_cerberus':
        fname = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/mousewg/refs/cerberus/vM25_gene_metadata.tsv'
    elif how in iso_hows and ver == 'vM25_cerberus':
        fname = '/Users/fairliereese/Documents/programming/mortazavi_lab/data/mousewg/refs/cerberus/vM25_transcript_metadata.tsv'


    df = pd.read_csv(fname, sep='\t')

    # pdb.set_trace()

    if how == 'gene':
        id_col = 'gid'
    elif how == 'iso':
        id_col = 'tid'
    elif how == 'ic':
        iso_col = 'tid'
        id_col = 'ic'
    elif how == 'tss':
        iso_col = 'tid'
        id_col = 'tss'
    elif how == 'tes':
        iso_col = 'tid'
        id_col = 'tes'

    # if using cerberus features, drop duplicate  entries that come
    # from the different transcripts using the same features
    # also ignore the transcript length column
    if how in ['ic', 'tss', 'tes']:
        df = add_feat(df, kind=id_col, col=iso_col)
        df.drop([iso_col, 't_len'], axis=1, inplace=True)

        # double check
        n_feats = len(df[id_col].unique().tolist())
        n_drop_dupes = len(df.drop_duplicates().index)
        if n_feats != n_drop_dupes:
            print('Warning: number of unique {}s does not match length of deduped info table'.format(how))
        df.drop_duplicates(inplace=True)

    if subset == 'polya':
        polya_cats = ['protein_coding', 'lncRNA', 'pseudogene']
        df = df.loc[df.biotype_category.isin(polya_cats)]
    elif subset == 'protein_coding':
        df = df.loc[df.biotype_category == 'protein_coding']
    elif subset == 'pseudogene':
        df = df.loc[df.biotype_category == 'pseudogene']
    elif subset == 'tf':
        df = df.loc[df.tf == True]

    biotype_counts = df[[id_col, 'biotype']].groupby('biotype').count()
    biotype_counts.reset_index(inplace=True)
    biotype_counts.rename({id_col: 'gencode_counts'}, axis=1, inplace=True)

    biotype_cat_counts = df[[id_col, 'biotype_category']].groupby('biotype_category').count()
    biotype_cat_counts.reset_index(inplace=True)
    biotype_cat_counts.rename({id_col: 'gencode_counts'}, axis=1, inplace=True)

    # add stable gid if requested
    if add_stable_gid:
        df['gid_stable'] = cerberus.get_stable_gid(df, 'gid')

    return df, biotype_counts, biotype_cat_counts

def get_tpm_table(df,
                    how='gene',
                    groupby='library',
                    nov=None,
                    min_tpm=None,
                    sample='all',
                    gene_subset=None,
                    save=False):
    """
    Parameters:
        df (pandas DataFrame): TALON abundance table
        how (str): Choose from 'gene' or 'iso'
        groupby (str): Choose from 'library' or 'sample'. Sample will avg.
        nov (list of str): List of accepted novelty types (w/ how='iso')
        min_tpm (float): Keep only genes / isos that have at least one
            TPM >= the value across the libraries
        gene_subset (str): Choose from 'polya' or None
        save (bool): Whether or not to save the output matrix

    Returns:
        df (pandas DataFrame): TPMs for gene or isoforms in the requested
            samples above the input detection threshold.
        ids (list of str): List of str indexing the table
    """
    print('Calculating {} TPM values'.format(how))

    dataset_cols = get_dataset_cols(df, sample=sample)
    df = rm_sirv_ercc(df)

    # merge with information about the gene
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'biotype_category', 'tf']]
    df = df.merge(gene_df, how='left', left_on='annot_gene_id', right_on='gid')

    # get indices that we'll need to subset on
    if how == 'gene':
        id_col = 'annot_gene_id'
        nov_col = 'gene_novelty'
        nov = ['Known']
    elif how == 'iso':
        id_col = 'annot_transcript_id'
        nov_col = 'transcript_novelty'

    # filter on novelty
    if nov:
        print('Subsetting for novelty categories {}'.format(nov))
        nov_inds = df.loc[df[nov_col].isin(nov), id_col].tolist()
    else:
        nov_inds = df[id_col].tolist()

    # filter on gene subset
    if gene_subset:
        print('Subsetting for {} genes'.format(gene_subset))
        if gene_subset == 'polya':
            polya_cats = ['protein_coding', 'pseudogene', 'lncRNA']
            gene_inds = df.loc[df.biotype_category.isin(polya_cats), id_col].tolist()
        elif gene_subset == 'tf':
            gene_inds = df.loc[df.tf == True, id_col].tolist()
    else:
        gene_inds = df[id_col].tolist()

    # get intersection of both
    subset_inds = list(set(nov_inds)&set(gene_inds))

    # sum up counts across the same gene
    if how == 'gene':
        df = df[dataset_cols+[id_col]]
        df = df.groupby(id_col).sum().reset_index()

    # set index so that all values in df reflect
    # counts per transcript or gene
    df.set_index(id_col, inplace=True)

    # compute TPM
    tpm_cols = []
    for d in dataset_cols:
        tpm_col = '{}_tpm'.format(d)
        total_col = '{}_total'.format(d)
        df[total_col] = df[d].sum()
        df[tpm_col] = (df[d]*1000000)/df[total_col]
        tpm_cols.append(tpm_col)
    df = df[tpm_cols]

    # reformat column names
    df.columns = [c.rsplit('_', maxsplit=1)[0] for c in df.columns]

    # enforce tpm threshold
    if min_tpm:
        print('Enforcing minimum TPM')
        print('Total # {}s detected: {}'.format(how, len(df.index)))
        df = df.loc[(df >= min_tpm).any(axis=1)]
        print('# {}s >= {} tpm: {}'.format(how, min_tpm, len(df.index)))

    # subset if necessary
    if gene_subset or nov:
        print('Applying gene type and novelty subset')
        df = df.loc[df.index.isin(subset_inds)]

    # average over biosample
    if groupby == 'sample':
        print('Averaging over biosample')
        df = df.transpose()
        df.reset_index(inplace=True)

        # add biosample name (ie without rep / sex information)
        df['biosample'] = df['index'].str.rsplit('_', n=2, expand=True)[0]

        print('Found {} total samples'.format(len(df.biosample.unique().tolist())))

        df = df.groupby('biosample').mean()
        df = df.transpose()

    print('Number of {}s reported: {}'.format(how, len(df.index)))

    if save:
        fname = '{}_{}_tpm.tsv'.format(sample, how)
        df.to_csv(fname, sep='\t')

    ids = df.index.tolist()

    return df, ids

def get_gene_number(df, col, pref):
    """
    Add number of tss / tes / ic occurrence to t_df
    from get_ic_tss_tes, based on number of unique
    values w/i the gene

    Parameters:
        df (pandas DataFrame): t_df from get_ic_tss_tes
        col (str): column with unique tss/ic/tes id
        pref (str): prefix to give new column

    Returns:
        df (pandas DataFrame): t_df with added numbers
            for tss/ic/tes id w/i gene
    """
    new_col = '{}_gene_num'.format(pref)
    temp = df[['gid', col, 'annotation']].copy(deep=True)
    temp.drop_duplicates(subset=['gid', col], inplace=True)
    temp[new_col] = temp.sort_values(['gid', col, 'annotation'],
                                 ascending=[True, True, False])\
                                 .groupby(['gid'])\
                                 .cumcount() + 1
    temp.drop('annotation', axis=1, inplace=True)
    df = df.merge(temp, how='left', on=['gid', col])
    return df

def count_tss_ic_tes(df, subset=None):
    """
    Count up unique tsss, ics, and tess for
    a given subset of transcript ids

    Parameters:
        df (pandas DataFrame): t_df from get_ic_tss_tes
        subset (list of str): List of transcript ids

    Returns:
        counts (pandas DataFrame): df w/ an entry detailing
            how many tss, ics, tes there are for each gene
    """
    df = df.copy(deep=True)
    df = df.loc[df.tid.isin(subset)]

    # raw tss, ic, tes count
    cols = ['tss', 'intron_chain', 'tes']
    for i, col in enumerate(cols):
        if col in ['tss', 'tes']:
            id_col = '{}_cluster'.format(col)
        else:
            id_col = col
        temp = df[[id_col, 'gid']].groupby('gid').nunique()
        if i == 0:
            counts = temp
        else:
            counts = counts.merge(temp, left_index=True, right_index=True)

    # unique combinations of tss, ic, tes
    df['tss_ic_tes'] = df.tss_cluster.astype(str)+'_'+\
                       df.intron_chain.astype(str)+'_'+\
                       df.tes_cluster.astype(str)
    temp = df[['tss_ic_tes', 'gid']].groupby('gid').nunique()
    counts = counts.merge(temp, how='left', left_index=True, right_index=True)

    for col in counts.columns:
        if col == 'tss_cluster':
            counts.rename({col: 'tss'}, axis=1, inplace=True)
        elif col == 'tes_cluster':
            counts.rename({col: 'tes'}, axis=1, inplace=True)

    # compute splicing ratio
    counts['splicing_ratio'] = counts.intron_chain/((counts.tes+counts.tss)/2)

    return counts

def df_to_pyranges(ends, kind='tss'):

    # reformat column names if needed
    cols = ends.columns
    if 'Start' in cols and 'End' in cols and 'Chromosome' in cols:
        pass
    else:
        coord = '{}_coord'.format(kind)
        chrom = '{}_chrom'.format(kind)
        ends = ends[cols].copy(deep=True)
        ends.rename({coord: 'Start',
                     chrom: 'Chromosome'},
                     axis=1, inplace=True)
        ends['End'] = ends.Start

    # turn into a pyranges object
    cols = ['gid', 'gname',
            'Start', 'End',
            'Chromosome', kind]
    if kind == 'tss':
        cols.append('first_sd')
    ends = ends[cols]
    ends.drop_duplicates(inplace=True)
    ends = pr.PyRanges(df=ends)

    return ends

def cluster_ends(ends,
                 slack,
                 cluster_start=1,
                 kind='tss'):
    """
    Cluster TSSs / TESs.

    Parameters:
        ends (pandas DataFrame): Slice of dataframe from add_tss_ic_tes
        slack (int): Allowable distance for merging ends
        cluster_start (int): # to start numbering clusters from
        kind (str): 'tss' or 'tes'

    Returns:
        reg (pandas DataFrame): DF describing found regions
        clust (pandas DataFrame): DF describing which region
            each end observation corresponds to
    """

    ends = df_to_pyranges(ends, kind=kind)

    # get regions and region assignments
    cols = ['gid', 'gname']
    if kind == 'tss':
        cols.append('first_sd')

    # merge to get regions
    reg = ends.merge(strand=None, by=cols, slack=slack)
    reg = reg.as_df()
    reg['len'] = reg.End - reg.Start
    reg['Cluster'] = [i for i in range(cluster_start, len(reg.index)+cluster_start)]

    # cluster to get region assignment per end
    clust = ends.cluster(strand=None, by=cols, slack=slack)
    clust = clust.as_df()
    clust['Cluster_new'] = clust.Cluster+cluster_start-1
    clust.drop('Cluster', axis=1, inplace=True)
    clust.rename({'Cluster_new': 'Cluster'}, axis=1, inplace=True)

    return reg, clust

def add_tss_ic_tes(sg):
    """
    Adds the intron chain, tss, tes, and first splice donor
    of each transcript to the t_df object. Also adds coords
    of tss and tes

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with annotated and observed transcripts

    Returns
        df (pandas DataFrame): DF with start / end vertex / coord
            info, ic, and first splice donor vertex info
    """
    df = sg.t_df.copy(deep=True)

    # add intron chains
    paths = df.path.values.tolist()
    paths = [tuple(path[1:-1]) for path in paths]
    df['intron_chain'] = paths

    # add tss
    paths = df.loc_path.values.tolist()
    tsss = [path[0] for path in paths]
    df['tss'] = tsss

    # add tes
    paths = df.loc_path.values.tolist()
    tess = [path[-1] for path in paths]
    df['tes'] = tess

    # add first splice donor
    paths = df.loc_path.values.tolist()
    first_sds = [path[1] for path in paths]
    df['first_sd'] = first_sds

    # add coordinates
    cols = ['tss', 'tes']
    for c in cols:
        # first, add tss / tes coords
        df = df.merge(sg.loc_df[['vertex_id', 'chrom', 'coord']],
                  how='left', left_on=c, right_index=True)
        df.drop(['vertex_id'], axis=1, inplace=True)
        df.rename({'chrom': '{}_chrom'.format(c),
                  'coord': '{}_coord'.format(c)},
                  axis=1, inplace=True)

    return df

def get_ic_tss_tes(sg,
                   df,
                   min_tpm=1,
                   sample='all',
                   gene_subset='polya',
                   annot_slack=200,
                   novel_slack=100,
                   verbose=False):
    """
    Extract information about annotaed and observed tss, tes,
    and intron chain usage from a SwanGraph t_df.

    Parameters:
        sg (swan_vis SwanGraph): SwanGraph with both annotation
            and observed transcript data added
        df (pandas DataFrame): Filtered TALON abundance file
        min_tpm (int): Min TPM to consider a transcript detected
            and therefore to include ends in end calling
        gene_subset (str): Choose from 'polya', 'tf'
        annot_slack (int): Distance b/w which to merge annotated ends
        novel_slack (int): Distance b/w which to merge observed ends
        verbose (bool): Whether or not to print output

    Returns:
        all_df (pandas DataFrame): sg.t_df modified to include
            information about intron chain, tss, and tes
        regions (dict of pandas DataFrames): Indexed by
            'tss' and 'tes'. Bed regions for each end cluster
            as annotated in all_df
        counts (pandas DataFrame): DF of counts for intron
            chains, TSSs, TESs, and unique combinations of the three
            for annotated, observed, and both
    """

    all_df = add_tss_ic_tes(sg)

    # limit to those annotated or in list of detected tids that we allow
    # only ever allow known, nic, nnc
    _, inds = get_tpm_table(df,
                             how='iso',
                             min_tpm=min_tpm,
                             sample=sample,
                             gene_subset=gene_subset,
                             nov=['Known', 'NIC', 'NNC'])
    novel_tids = inds

    if type(novel_tids) == list:
        all_df = all_df.loc[(all_df.annotation == True)|(all_df.tid.isin(novel_tids))]

    end_types = ['tss', 'tes']
    end_regions = dict()
    for c in end_types:

        if verbose:
            print()

        #### annotated transcripts ####

        t_df = all_df.loc[all_df.annotation == True].copy(deep=True)
        if verbose:
            n = len(t_df.index)
            print('Finding {}s for {} annotated transcripts'.format(c, n))

        reg, clust = cluster_ends(t_df,
                                  slack=annot_slack,
                                  cluster_start=1,
                                  kind=c)
        reg['annotation'] = True
        reg['source'] = 'GENCODE'
        clust['annotation'] = True
        if verbose:
            n = len(reg.index)
            print('Found {} annotated {} clusters'.format(n,c))

        #### novel transcripts ####

        # assign ends from novel transcripts to annotated ends
        t_df = all_df.loc[(all_df.annotation == False)&(all_df.tid.isin(sg.adata.var.index.tolist()))]

        if verbose:
            n = len(t_df.index)
            print('Finding {}s for {} novel transcripts'.format(c, n))

        # case 1: ends from novel transcripts are w/i annotated regions
        ends = df_to_pyranges(t_df, kind=c)
        reg = pr.PyRanges(df=reg)
        ends = ends.join(reg, how='left',
                         slack=1,
                         strandedness=None, suffix='_annot')
        ends = ends.as_df()

        # limit to those w/ matching gid, first sd and add to cluster df
        if c == 'tss':
            inds = ends.loc[(ends.gid==ends.gid_annot)&(ends.first_sd==ends.first_sd_annot)].index.tolist()
        else:
            inds = ends.loc[ends.gid == ends.gid_annot].index.tolist()

        if verbose:
            n = len(inds)
            print('Found {} novel {}s that are already in the annotation'.format(n,c))
        clust = pd.concat([clust, ends.loc[inds]])

        # case 2: ends from novel transcripts need to be clustered
        # on their own

        # remove duplicates that arise from imperfect merging
        ends['in_region'] = False
        ends.loc[inds, 'in_region'] = True
        ends.sort_values(by='in_region', inplace=True, ascending=False)
        cols = ['gid', c]
        if c == 'tss':
            cols.append('first_sd')
        ends.drop_duplicates(subset=cols, keep='first', inplace=True)
        inds = ends.loc[ends.in_region == True].index.tolist()

        # subset based on ends that are unsupported by ref regions
        inds = list(set(ends.index.tolist())-set(inds))
        t_df = ends.loc[inds]
        if verbose:
            n = len(t_df.index)
            print('Finding {}s for {} novel ends'.format(c,n))
        n = clust.Cluster.max()+1
        nov_reg, nov_clust = cluster_ends(t_df,
                                          slack=novel_slack,
                                          cluster_start=n,
                                          kind=c)
        nov_reg['annotation'] = False
        nov_reg['source'] = 'obs'
        nov_clust['annotation'] = False
        if verbose:
            n = len(nov_reg.index)
            print('Found {} novel {} clusters'.format(n,c))

        # check how many novel clusters fall into already
        # annotated regions
        nov_reg = pr.PyRanges(df=nov_reg)
        temp = nov_reg.join(reg, how=None, strandedness=None, suffix='_annot')
        temp = temp.as_df()
        if verbose:
            if c == 'tss':
                temp = temp.loc[(temp.first_sd == temp.first_sd_annot)&(temp.gid == temp.gid_annot)]
            else:
                temp = temp.loc[temp.gid == temp.gid_annot]
            cols = ['gid']
            if c == 'tss':
                cols.append('first_sd')
            temp = temp.drop_duplicates(subset=cols)
            n = len(temp.index)
            print('{} new {} regions overlap annotated regions'.format(n,c))

        # finally, add novel regions to clust and reg dfs
        reg = reg.as_df()
        nov_reg = nov_reg.as_df()
        clust = pd.concat([clust, nov_clust])
        reg = pd.concat([reg, nov_reg])

        # some final formatting for these dfs
        cols = ['gid', 'gname', c,
                'Cluster', 'annotation']
        if c == 'tss':
            cols.append('first_sd')
        clust = clust[cols]
        clust.rename({'Cluster': '{}_cluster'.format(c),
                      'annotation': '{}_annotation'.format(c)},
                      axis=1, inplace=True)
        clust.drop_duplicates(inplace=True)
        end_regions[c] = reg

        # add cluster ids back into the original df
        cols = ['gid', 'gname', c]
        if c == 'tss':
            cols.append('first_sd')
        all_df = all_df.merge(clust, how='left', on=cols)

    # counts for all, annotated, and observed go into the same df
    # with a different source
    counts = pd.DataFrame()

    # annotated counts
    tids = all_df.loc[all_df.novelty == 'Known'].tid.tolist()
    annot_counts = count_tss_ic_tes(all_df, subset=tids)
    annot_counts['source'] = 'GENCODE'
    counts = pd.concat([counts, annot_counts])

    # annotated + novel counts
    tids = list(set(tids)|set(novel_tids))
    all_counts = count_tss_ic_tes(all_df, subset=tids)
    all_counts['source'] = 'all'
    counts = pd.concat([counts, all_counts])

    # observed counts
    tids = list(set(novel_tids)&set(sg.adata.var.index.tolist()))
    obs_counts = count_tss_ic_tes(all_df, subset=tids)
    obs_counts['source'] = 'obs'
    counts = pd.concat([counts, obs_counts])

    # get gene info and add
    gene_df, _, _ = get_gtf_info(how='gene')
    gene_df = gene_df[['gid', 'gname', 'biotype',
                  'biotype_category', 'tf']]
    gene_df.drop_duplicates(inplace=True)
    counts = counts.merge(gene_df, how='left', left_index=True, right_on='gid')

    t_df = all_df.copy(deep=True)

    # add tripletized transcript name to t_df
    t_df = get_gene_number(t_df, 'tss_cluster', 'tss')
    t_df = get_gene_number(t_df, 'tes_cluster', 'tes')
    t_df = get_gene_number(t_df, 'intron_chain', 'intron_chain')
    t_df['ttrip'] = t_df.gname +' ('+\
                t_df.tss_gene_num.astype('str')+','+\
                t_df.intron_chain_gene_num.astype('str')+','+\
                t_df.tes_gene_num.astype('str')+')'

    return t_df, end_regions, counts
