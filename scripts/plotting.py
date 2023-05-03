import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
# from matplotlib_venn import venn2
from matplotlib.colors import ListedColormap
import matplotlib as mpl
# import upsetplot
from scipy import stats
import os
import matplotlib.colors as mc
# import ternary
from sklearn import preprocessing
import pylab as pl
import matplotlib.ticker as tck
from collections import defaultdict

# from .utils import *

def get_lr_bulk_sample_colors():
    c_dict, order = get_tissue_age_colors()

    # c2c12
    c_dict['c2c12_myoblast'] = '#ca79a7'
    c_dict['c2c12_myotube'] = '#009c73'
    
    # forelimb
    c_dict['forelimb_e11'] = '#99ebec'
    c_dict['forelimb_e13'] = '#01ced0'
    
    # adrenal, hc, ctx
    for t in ['adrenal', 'hippocampus', 'cortex']:
        c_dict[t] = get_tissue_colors()[0][t]   
        
    return c_dict, None

def get_sector_colors():
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    simple = '#e5ecf6'
    c_dict = {'tss': tss,
              'splicing': splicing,
              'tes': tes,
              'simple': simple}
    order = ['tss', 'splicing', 'tes', 'simple']
    return c_dict, order

def get_tissue_age_colors():
    c_dict, order = get_tissue_colors()

    # get different color for each age / tissue
    new_c_dict = {}
    min_color = '#FFFFFF'
    ages = ['4d', '10d', '14d', '25d', '36d', '2mo', '18-20mo']
    order = []
    for tissue, max_color in c_dict.items():
        cmap = mpl.colors.LinearSegmentedColormap.from_list(tissue, [min_color, max_color], N=8)
        for i, age in enumerate(ages):
            key = '{}_{}'.format(tissue, age)
            new_c_dict[key] = mpl.colors.to_hex(cmap(i+1))
            order.append(key)

    return new_c_dict, order

def get_celltype_colors(dataset, res):
    """
    Get celltype colors for different tissues at different resolutions
    
    Parameters:
        dataset (str): {'adrenal', 'hippocampus', 'cortex'
                        'gastrocnemius', 'heart'}
        res (str): {'celltype', 'sub_celltype', 'gen_celltype'}
    """
    
    d = os.path.dirname(__file__)
    fname = '{}/../refs/enc4_mouse_snrna_celltypes_c2c12.csv'.format(d)
    df = pd.read_csv(fname)
    
    tissue = dataset
    tissue = tissue.capitalize()
    if res == 'celltype':
        df_cat = 'celltypes'
        df_color = 'celltype_color'
    elif res == 'sub_celltype':
        df_cat = 'subtypes'
        df_color = 'subtype_color'
    elif res == 'gen_celltype':
        df_cat = res
        df_color = 'gen_celltype_color'

    temp = df.loc[df.tissue == tissue]
    temp = temp[[df_cat, df_color]]
    temp.drop_duplicates(inplace=True)
    temp.set_index(df_cat, inplace=True)
    temp.head()
    c_dict = temp.to_dict()[df_color]
    return c_dict

# def get_celltype_colors(dataset='adrenal', cats=None, rm_cats=None):
#     """
#     Return a dictionary of celltype : color as well as order of
#     celltypes for the given dataset.

#     Parameters:
#         dataset (str): Choose from 'adrenal', 'hippocampus'
#         cats (list of str): List of celltypes to include
#         rm_cats (list of str): List of celltypes to exclude
#     """
#     if dataset == 'adrenal':
#         order =  ['Medulla_NE','Medulla_EPI','Sox10+',
#                   'Stromal','Adipocytes','Hepatocyte',
#                   'Smooth_muscle','Macrophage','Endothelial',
#                   'Cortex/Endothelial','Cortex_ZF','Cortex_ZG',
#                   'Cortex_cycling','X_zone','Y_zone',
#                   'Capsule','Other']

#         c_dict = {'Medulla_NE':'#97578a',
#                     'Medulla_EPI':'#339470',
#                     'Sox10+':'#753b74',
#                     'Stromal':'#e2969b',
#                     'Adipocytes':'#e25e2c',
#                     'Hepatocyte':'#ad7797',
#                     'Smooth_muscle':'#86b84d',
#                     'Macrophage':'#da5774',
#                     'Endothelial':'#a63b4c',
#                     'Cortex/Endothelial':'#f99d26',
#                     'Cortex_ZF':'#f0c130',
#                     'Cortex_ZG':'#b2373a',
#                     'Cortex_cycling':'#3f4075',
#                     'X_zone':'#fade7c',
#                     'Y_zone':'#e47381',
#                     'Capsule':'#236d88',
#                     'Other':'grey'}
#     return c_dict, order

def get_tissue_colors(cats=None, rgb=False):
    d = os.path.dirname(__file__)
    fname = '{}/../refs/tissue_colors.tsv'.format(d)
    df = pd.read_csv(fname, sep='\t')
    c_dict = {}
    for ind, entry in df.iterrows():
        c_dict[entry.tissue] = entry.color
    order = ['adrenal', 'hippocampus', 'cortex', 'gastroc', 'heart']

    if cats:
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]

    if rgb:
        for key, item in c_dict.items():
            item = item[1:]
            r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
            c_dict[key] = (r,g,b)

    return c_dict, order

def get_age_colors():
    """
    Return a dictionary mapping timepoints to colors
    """

    c_dict = {'4d': 'thistle',
             '10d': 'plum',
             '14d': 'violet',
             '25d': 'mediumorchid',
             '36d': 'darkviolet',
             '2mo': 'purple',
             '18-20mo': 'indigo'}
    order = ['4d', '10d', '14d', '25d',
             '36d', '2mo', '18-20mo']
    return c_dict, order

def get_sex_colors():
    """
    Return a dictionary mapping sexes to colors
    """
    c_dict = {'m': '#e47381', 'f': '#fade7c'}
    order = ['f', 'm']
    return c_dict, order

# def get_tissue_age_colors(cats=None):
#     d = os.path.dirname(__file__)
#     fname = '{}/../refs/tissue_age_colors.tsv'.format(d)
#     df = pd.read_csv(fname, sep='\t')

#     order = df.tissue_age.tolist()
#     c_dict = None

#     return c_dict, order

def get_talon_nov_colors(cats=None):
    c_dict = {'Known': '#009E73',
              'ISM': '#0072B2',
              'NIC': '#D55E00',
              'NNC': '#E69F00',
              'Antisense': '#000000',
              'Intergenic': '#CC79A7',
              'Genomic': '#F0E442'}
    order = ['Known', 'ISM', 'NIC', 'NNC', 'Antisense', 'Intergenic', 'Genomic']

    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]
    return c_dict, order

def set_metadata_colors(adata, obs_col, cmap):
		"""
		Set plotting colors for datasets based on a column in the metadata
		table.
		Parameters:
			obs_col (str): Name of metadata column to set colors for
			cmap (dict): Dictionary of metadata value : color (hex code with #
				character or named matplotlib color)
		"""
		# pdb.set_trace()

		# map values in order specific to
		adata.obs[obs_col] = adata.obs[obs_col].astype('category')
		obs_order = list(adata.obs_names)
		sample_order = adata.obs[obs_col].cat.categories.tolist()
		sample_colors = [cmap[s] for s in sample_order]
		adata.uns['{}_colors'.format(obs_col)] = sample_colors

		# if colors are named, get hex code
		for key, item in cmap.items():
			if '#' not in item:
				cmap[key] = mc.cnames[item]

# 		# also store rgb values in dict for use with gen_report
# 		for key, item in cmap.items():
# 			item = item[1:]
# 			r,g,b = tuple(int(item[i:i+2], 16) for i in (0, 2, 4))
# 			cmap[key] = (r,g,b)
# 		adata.uns['{}_dict'.format(obs_col)] = cmap

def plot_short_long_det(df, opref, \
                    xlim, ylim, how='gene'):

    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'ill_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'ill_umi_count'

#     ax = sns.jointplot(data=df, x=c1, y=c2,
#                      xlim=(0,xlim), ylim=(0,ylim),
#                      joint_kws={'data':df, 's':40, 'alpha':1})
    ax = sns.jointplot(data=df, x=c1, y=c2,
                     xlim=(0,xlim), ylim=(0,ylim),
                     joint_kws={'data':df, 's':40, 'alpha':1})
    ax = ax.ax_joint

#     # plot regression lines and equation of regression lines
#     # https://stackoverflow.com/questions/48145924/different-colors-for-points-and-line-in-seaborn-regplot/68135585#68135585
#     # https://stackoverflow.com/questions/45902739/seaborn-annotate-the-linear-regression-equation
#     # https://stackoverflow.com/questions/62705904/add-entry-to-matplotlib-legend-without-plotting-an-object
#     lines = []
#     labels = []
#     for s in df['sample'].unique().tolist():
#         temp = df.loc[df['sample'] == s]
#         color = c_dict[s]
#         line_color = adjust_lightness(color, 0.5)

#         # get coeffs of linear fit
#         slope, intercept, r_value, p_value, std_err = stats.linregress(temp[c1],temp[c2])
#         lines += [mpl.lines.Line2D([0], [0], color=line_color)]
#         labels += ['m={0:.1f}'.format(slope)]

#         print('Slope of {} correlation: {}'.format(s, slope))

#         sns.regplot(data=temp, x=c1, y=c2,
#                     scatter=False, ax=ax, color=color)
#         sns.regplot(data=temp, x=c1, y=c2,
#             scatter=False, ax=ax, color=color, ci=0,
#             line_kws={'color':line_color,
#                       'linestyle':'-',
#                       'label':"m={0:.1f}".format(slope)})

    ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')


def plot_reads_per_bc(df, title, oprefix):
    """
    Parameters:
        df (pandas DataFrame): DataFrame of read_annot file
        title (str): Title of plot
        oprefix (str): Output file prefix
    """

    temp = get_reads_per_bc(df)

    sns.set_context('paper', font_scale=1.5)
    counts = temp['counts'].tolist()
    plt.plot(range(len(counts)),
            counts,
            color='lightgray',
            linewidth=2)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_xlabel('Ranked cells by # reads (logscale)')
    ax.set_ylabel('# reads (logscale)')
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    if title == 'Post-TALON':
        title = 'post_talon'

    plt.tight_layout()

    fname = '{}_{}_umis_v_barcodes.png'.format(oprefix, title)
    plt.savefig(fname)

def plot_blastp_scores(df, opref, hue=None,
                      order=None, c_dict=None):
    sns.set_context('paper', font_scale=2)

    if not hue:
        ax = sns.displot(data=df, x='blastp_score', kind='kde', linewidth=3)
    else:
        if c_dict and order:
            ax = sns.displot(data=df, x='blastp_score', kind='kde', linewidth=3, hue=hue, palette=c_dict, hue_order=order, common_norm=False)
        else:
            ax = sns.displot(data=df, x='blastp_score', kind='kde', linewidth=3, hue=hue, common_norm=False)
    if hue:
        fname = '{}_{}_blastp_score.pdf'.format(opref, hue)
    else:
        fname = '{}_blastp_score.pdf'.format(opref)
    plt.savefig(fname, dpi=300, bbox_inches='tight')

def plot_read_len_kde(df, opref, hue=None, common_norm=False):
    sns.set_context('paper', font_scale=2)

#     ax = sns.displot(data=df, x='read_length', hue=hue,
#                  palette=c_dict, kind='kde', hue_order=order, linewidth=3,
#                  common_norm=common_norm)
    if not hue:
        ax = sns.displot(data=df, x='read_length', kind='kde',
                         linewidth=3, common_norm=common_norm)
    else:
        ax = sns.displot(data=df, x='read_length', kind='kde',
                         linewidth=3, common_norm=common_norm, hue=hue)
#     ax.set(xlabel='Read length', ylabel='KDE of reads',
#           title='Length distribution of Reads', xlim=(0,7500),
#           xticks=[0, 2500, 5000, 7500])
    plt.savefig('{}_read_length_kde.pdf'.format(opref), dpi=300, bbox_inches='tight')

def plot_cluster_proportions(adata,
                             bar_key,
                             color_key,
                             scale=False,
                             opref='figures/'):
    """
    Parameters
        scale (bool): Whether or not to scale props.
            to relative number of cells
    """

    # raw proportions
    if not scale:
        adata_tmp = adata.copy()
        sizes = adata_tmp.obs.groupby([color_key, bar_key]).size()
        props = sizes.groupby(level=1).apply(lambda x: 100 * x / x.sum()).reset_index()
        props = props.pivot(columns=bar_key, index=color_key).T
        props.index = props.index.droplevel(0)
        props.fillna(0, inplace=True)

    # scaled
    else:
        temp2 = adata.obs.copy(deep=True)
        temp2.index.name = 'index'
        temp2.reset_index(inplace=True)
        temp2 = temp2[['index', color_key]].groupby(color_key).count().reset_index()
        temp2.rename({'index': 'color_counts'}, axis=1, inplace=True)

        temp = adata.obs.copy(deep=True)
        temp.index.name = 'index'
        temp.reset_index(inplace=True)

        temp = temp[['index', color_key, bar_key]].groupby([color_key, bar_key]).count().reset_index()
        temp.rename({'index': 'counts'}, axis=1, inplace=True)
        temp = temp.merge(temp2, how='left', on='experiment')
        temp['prop_color'] = temp.counts/temp.color_counts

        temp3 = temp[[bar_key, 'prop_color']].groupby(bar_key).sum().reset_index()
        temp3.rename({'prop_color': 'scale_factor'}, axis=1, inplace=True)

        temp = temp.merge(temp3, how='left', on=bar_key)
        temp['scaled_prop'] = temp.prop_color/temp.scale_factor

        temp = temp[[bar_key, color_key, 'scaled_prop']]
        temp = temp.pivot(index=bar_key, columns=color_key, values='scaled_prop')
        props = temp

    fig, ax = plt.subplots(dpi=300)
    fig.patch.set_facecolor("white")
    fig.set_size_inches(15, 5)

    cmap = None
    cluster_palette = '{}_colors'.format(color_key)
    if cluster_palette in adata.uns.keys():
        cluster_palette = adata.uns[cluster_palette]
        cmap = sns.palettes.blend_palette(
            cluster_palette,
            n_colors=len(cluster_palette),
            as_cmap=True)

    props.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        legend=None,
        colormap=cmap
    )

    ax.legend(bbox_to_anchor=(1.1, 2), title=color_key)
    sns.despine(fig, ax)
    ax.tick_params(axis="x", rotation=90)
    ax.set_xlabel(props.index.name.capitalize())
    if not scale:
        ax.set_ylabel("Percent")
    else:
        ax.set_ylabel('Scaled proportion')
    ax.grid(False)

    fname = opref+'{}_by_{}_prop.png'.format(bar_key, color_key)
    plt.savefig(fname, bbox_inches='tight')
    
    return props

def plot_depth_by_tech(adata, how, opref,
                       hue=None, xlim=None, ylim=None):
    df = adata.obs
    sns.set_context('paper', font_scale=2)

    if how == 'gene':
        c1 = 'n_genes'
        c2 = 'sr_gene_count'
    elif how == 'read':
        c1 = 'n_counts'
        c2 = 'sr_umi_count'

    if xlim:
        xlim = (0, xlim)
    if ylim:
        ylim = (0, ylim)
    ax = sns.jointplot(data=df, x=c1, y=c2,
                       xlim=xlim, ylim=ylim,
                       joint_kws={'data':df, 's':40, 'alpha':1},
                       hue=hue)
    ax = ax.ax_joint

    # ax.legend(title='')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    # ax.get_legend().remove()

    if how == 'gene':
        _ = ax.set(xlabel='# genes/cell in PacBio', ylabel='# genes/cell detected in Illumina')
        plt.savefig('{}_genes_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')
    elif how == 'read':
        _ = ax.set(xlabel='# reads/cell in PacBio', ylabel='# UMIs/cell in Illumina')
    plt.savefig('{}_reads_detected_pb_v_illumina.pdf'.format(opref), dpi=300, bbox_inches='tight')

def add_perc(ax, data, feature):
    total = data[feature].sum()
    ylim = ax.get_ylim()[1]
    for p in ax.patches:
        percentage = '{:.1f}%'.format(100 * p.get_height()/total)
        x = p.get_x() + p.get_width() / 2 - 0.45
        y = p.get_y() + p.get_height() + ylim*0.00625
        ax.annotate(percentage, (x, y), size = 12)

def plot_read_novelty(df, opref, c_dict, order,
                      ylim=None, title=None,
                      datasets='all'):
    sns.set_context("paper", font_scale=1.6)

    temp = df.copy(deep=True)

    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]

    # count number of reads per cat
    temp = temp[['transcript_novelty', 'read_name']].groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'read_name':'counts'}, axis=1, inplace=True)
    print(temp)

    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Reads')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

    # add title
    if not title:
        g.fig.suptitle('Reads per novelty category')
    else:
        g.fig.suptitle('{} reads per novelty category'.format(title))

    # save figure
    fname = '{}_read_novelty'.format(opref)
    g.savefig(fname+'.pdf', dpi=300)

def plot_transcript_novelty(df, oprefix, c_dict, order, \
                            ylim=None, title=None,
                            whitelist=None, datasets='all', save_type='pdf'):
    sns.set_context('paper', font_scale=1.6)

    temp = df.copy(deep=True)

    # remove transcripts that are not on whitelist
    if whitelist:
        temp = temp.loc[temp.transcript_ID.isin(whitelist)]

    # filter on datasets
    if datasets != 'all':
        temp = temp.loc[temp.dataset.isin(datasets)]

    # count number of isoforms per cat
    temp = temp[['transcript_ID', 'transcript_novelty', 'read_name']].groupby(['transcript_ID', 'transcript_novelty']).count()
    temp.reset_index(inplace=True)
    temp.drop('read_name', axis=1, inplace=True)
    temp = temp.groupby('transcript_novelty').count()
    temp.reset_index(inplace=True)
    temp.rename({'transcript_ID': 'counts'}, axis=1, inplace=True)
    print(temp)

    # actual plotting
    g = sns.catplot(data=temp, x='transcript_novelty',
                y='counts', kind='bar',
                palette=c_dict, order=order)
    [plt.setp(ax.get_xticklabels(), rotation=90) for ax in g.axes.flat]
    g.set_ylabels('Isoforms')
    g.set_xlabels('Transcript novelty')

    # add percentage labels
    ax = g.axes[0,0]
    add_perc(ax, temp, 'counts')

    if ylim:
        g.set(ylim=(0,ylim))

    # add title
    if not title:
        g.fig.suptitle('Transcript models per novelty category')
    else:
        g.fig.suptitle('{} transcript models per novelty category'.format(title))

    # save figure
    fname = '{}_isoform_novelty'.format(oprefix)
    if save_type == 'png':
        g.savefig(fname+'.png', dpi=300)
    elif save_type == 'pdf':
        g.savefig(fname+'.pdf', dpi=300)

    plt.show()
    plt.clf()

def zip_pts(df, c):
    return zip(df[c['a']], df[c['b']], df[c['c']])

def max_pts(df, c):
    return max(df[c['a']].max(), df[c['b']].max(), df[c['c']].max())

def density_dorito(counts,
                   c,
                   scale=20,
                   cmap='viridis',
                   vmax=None,
                   log=False,
                   pad=0.15):
    """
    Plot the density of a dataset on a ternary plot
    From here: https://github.com/marcharper/python-ternary/issues/81
    
    Parameters:
        counts 
        c
        scale 
        cmap 
        
    Returns: 
        fig
        tax
        counts (pandas DataFrame): Counts, scaled by factor used
    """
    hm_dict = defaultdict(int)
    for i in range(0, scale+1):
        for j in range(0, scale+1):
            for k in range(0, scale+1):
                if i+j+k == scale:
                    # print(i,j,k)
                    # i = 1
                    # j = 1
                    # k = 1
                    temp = counts.copy(deep=True)
                    if i != scale:
                        temp = temp.loc[(temp.tss_ratio*scale>=i)&(temp.tss_ratio*scale<i+1)]
                        # print(temp)
                    else:
                        temp = temp.loc[(temp.tss_ratio*scale>=i)&(temp.tss_ratio*scale<=i+1)]
                    if j != scale:
                        temp = temp.loc[(temp.top_ratio*scale>=j)&(temp.top_ratio*scale<j+1)]
                        # print(temp)
                    else:
                        temp = temp.loc[(temp.top_ratio*scale>=j)&(temp.top_ratio*scale<=j+1)]
                    # print(i)
                    # print(j)
                    # print(temp.head())
                    n = len(temp.index)
                    # print(n)
                    hm_dict[i,j] += n
            
    # log values if necessary
    if log:
        for key, item in hm_dict.items():
            hm_dict[key] = np.log2(item+1)
    
    # double checking stuff
    df = pd.DataFrame.from_dict(hm_dict, orient='index')
    df['i'] = [b[0] for b in df.index.tolist()]
    df['j'] = [b[1] for b in df.index.tolist()]
    # df['k'] = [b[2] for b in df.index.tolist()]
    df.rename({0:'val'}, axis=1, inplace=True)
    # print(df.loc[df.val >= 14])
    
    
        
    figure, tax = ternary.figure(scale=scale, permutation='210')
    # tax.heatmap(hm_dict, colorbar=False, style='t', vmax=vmax)
    tax.heatmap(hm_dict, colorbar=False, style='t', adj_vlims=True, cmap=cmap)
    # tax.heatmap(interp_dict, colorbar=False)
    
    # scale according to chosen resolution
    for key in c.keys():
        counts[c[key]] = counts[c[key]]*scale
        
    # colorbar - hacked together by broken ternary code 
    ax = tax.get_axes()
    flat = []
    for key, item in hm_dict.items():
        flat.append(item)
    min_val = min(flat)
    max_val = max(flat)
    
    if vmax: 
        max_val = vmax
        
    # print(min_val)
    # print(max_val)
    norm = plt.Normalize(vmin=min_val, vmax=max_val)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    
    def exp_format(x,pos):
            x = int(x)
            return r'$2^{{{}}}$'.format(x)
    
    if not log:
        cb = plt.colorbar(sm, ax=ax, pad=pad)
    else:
        cb = plt.colorbar(sm, ax=ax, pad=pad, 
                          format=tck.FuncFormatter(exp_format))
    
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(16)
    if not log:
        cb.ax.set_yticklabels([])
        cb.ax.set_yticks([])

    cb.set_label('Density', size=16)
    
    return figure, tax, counts

def rm_color_cats(c_dict, order, cats):
    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats] 
    return c_dict, order

def jitter_dorito(counts, c, scale):
    """
    Parameters:
        counts
        c
        scale

    Returns
        counts
        c
    """

    # figure out how much to jitter by
    sigma = (1/250)*scale
    for d in c.keys():
        d_jitter = '{}_jitter'.format(d)
        counts[d_jitter] = counts[c[d]].apply(lambda x: np.random.normal(x, sigma))
        c[d] = d_jitter

    return counts, c

def scatter_dorito(counts,
                   c,
                   hue,
                   size,
                   log_size,
                   cmap,
                   mmap,
                   alpha,
                   density,
                   legend,
                   figure,
                   tax):
    """
    Parameters
        counts (pandas DataFrame): subset the thing
        c (dict of str): Dictionary of column names to plot as a, b, c
            indexed by 'a', 'b', 'c'
    """

    def scale_col(points, counts, col, log=False, how='color'):
            if log:
                log_col = '{}_log'.format(col)
                counts[log_col] = np.log10(counts[col])
                col = log_col
            vals = counts[[col]]
            max_val = vals[col].max()
            min_val = vals[col].min()
            min_max_scaler = preprocessing.MinMaxScaler(feature_range=(10, 300))
            vals = min_max_scaler.fit_transform(vals)
            max_scaled = max(vals)
            min_scaled = min(vals)

            # replace nans w/ 100
            vals = [100 if np.isnan(v) else v for v in vals]

            return vals, min_val, max_val, min_scaled, max_scaled

    # defaults
    points = [(x[0], x[1], x[2]) for x in zip_pts(counts, c)]
    labels = ['' for i in range(len(points))]
    hue_type = None
    figsize = (10,10)
    colors = '#e78424'
    if len(points) < 60:
        sizes = [100 for i in range(len(points))]
    else:
        sizes =  [20 for i in range(len(points))]
    markers = 'o'
    vmin = 0
    vmax = 1
    plotted = False

    # get color
    if hue:

        # categorical
        if counts[hue].dtype.name == 'object':
            hue_type = 'cat'
            colors = counts[hue].map(cmap).tolist()
            labels = counts[hue].tolist()

        # continuous
        else:
            hue_type = 'cont'
            colors, abs_min, abs_max, vmin, vmax = scale_col(points, counts, hue)

    # get sizes
    if size:
        sizes,_,_,_,_ = scale_col(points, counts, size, log_size)
        print(sizes[:5])

    # marker style
    if mmap:
        markers = [mmap[val] if val in mmap.keys() else 'o' for val in counts[hue].unique()]

    # figure size handling
    if hue_type == 'cat' and density: figsize = (13,10)
    elif hue_type == 'cat' and not density: figsize = (10,10)
    elif hue_type == 'cont' and density: figsize = (16,10)
    elif hue_type == 'cont' and not density: figsize = (13,10)
    elif density: figsize = (13,10)
    figure.set_size_inches(figsize[0], figsize[1])

    # actual scatter call
    if hue_type == 'cat':
        for point, color, size, label, marker in zip(points, colors, sizes, labels, markers):
            tax.scatter([point], vmin=vmin, vmax=vmax,
                    s=size, c=color, cmap=cmap,
                    marker=marker,label=label,
                    alpha=alpha, zorder=3)
    else:
        tax.scatter(points, vmin=vmin, vmax=vmax,
                    s=sizes, c=colors, cmap=cmap, marker=markers,
                    alpha=alpha, zorder=3)

    # legend handling
    if hue_type == 'cat' and legend:
        if density: x = 1.6
        else: x = 1.4
        tax.legend(bbox_to_anchor=(x, 1.05),
                   loc='upper right', prop={'size': 14})

        # fix marker size
        ax = tax.get_axes()
        lgnd = ax.get_legend()
        for handle in lgnd.legendHandles:
            handle._sizes = [100]

    # colorbar handling
    if hue_type == 'cont':
        ax = tax.get_axes()
        norm = plt.Normalize(vmin=abs_min, vmax=abs_max)
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        cb = plt.colorbar(sm, ax=ax, pad=0.1)
        for t in cb.ax.get_yticklabels():
             t.set_fontsize(16)
        if hue == 'tss' or hue == 'tes':
            cb.set_label('# {}s'.format(hue.upper()), size=16)
        elif hue == 'intron_chain':
            cb.set_label('# {}s'.format(hue), size=16)

    return figure, tax

def line_dorito(alpha, beta, gamma,
                scale, tax, figure):
    c_dict, _ = get_sector_colors()

    # scale
    alpha = alpha*scale
    beta = beta*scale
    gamma = gamma*scale

    linewidth = 3

    # splicing line
    tax.horizontal_line(beta, linewidth=linewidth,
                        color=c_dict['splicing'],
                        linestyle='--')

    # tss
    tax.right_parallel_line(alpha, linewidth=linewidth,
                           color=c_dict['tss'],
                           linestyle='--')

    # tes
    tax.left_parallel_line(gamma, linewidth=linewidth,
                           color=c_dict['tes'],
                           linestyle='--')

def plot_dorito(counts,
                top='splicing_ratio',
                subset=None,
                gene=None,
                hue=None,
                cmap='magma',
                mmap=None,
                density=False,
                density_scale=1,
                density_cmap='viridis',
                density_vmax=None,
                sectors=False,
                sect_alpha=0.5,
                sect_beta=0.5,
                sect_gamma=0.5,
                log_density=False,
                scatter=True,
                size=None,
                legend=True,
                log_size=False,
                jitter=False,
                alpha=1,
                scale=True,
                title=None,
                opref='figures/'):
    """
    Plot a dorito from counts with the given subset in a given
    color

    Parameters:
        counts (pandas DataFrame): DF of the counts per gene
            of ic, tss, tes from get_ic_tss_tes or
            compute_triplets (or both!!!)
        top (str): Column name to plot as apex of dorito.
            Choose from 'ic' or 'splicing_ratio'
        subset (dict of lists): List mapping counts column names
            to values in said columns to include in the data
        hue (str): Column from counts to color by
        cmap (str or dict of str): Either a dictionary mapping
            categorical column values from hue or a valid
            matplotlib continuous named color palette
        mmap (str or dict of str): Dictionary mapping categorical
            column values from hue to marker styles
        scale (bool): Whether to scale values b/w 1 and 0.
        alpha (float): Alpha value of points
        title (str): Title to give plot
        opref (str): Output file prefix to save fig
    """

    #### subset dataset and transform numbers as needed ####
    temp = counts.copy(deep=True)

    # if we have a gene name, limit to those entries
    if gene:
        temp = temp.loc[temp.gname == gene]

    # if we have a list of allowed sources, limit to those entries
    if subset:
        for col, val in subset.items():
            if type(val) != list:
                val = [val]
            temp = temp.loc[temp[col].isin(val)]

    # scale and assign which columns to use
    c = dict()
    if scale:
        if top == 'splicing_ratio':
            temp['total'] = temp.n_tss+temp.n_tes+temp.splicing_ratio
        elif top == 'intron_chain':
            temp['total'] = temp.n_tss+temp.n_tes+temp.intron_chain
        temp['tss_ratio'] = temp.n_tss/temp.total
        temp['tes_ratio'] = temp.n_tes/temp.total
        temp['top_ratio'] = temp[top]/temp.total

        c['a'] = 'tss_ratio'
        c['b'] = 'top_ratio'
        c['c'] = 'tes_ratio'
    else:
        c['a'] = 'tss'
        c['b'] = top
        c['c'] = 'tes'

    if scale == True:
        scale = 1
        mult = 0.2
    else:
        scale = max_pts(temp, c)

    # density
    if density:
        if hue:
            if counts[hue].dtype.name == 'object':
                pad = 0.1
            else:
                pad = 0.0
        else:
            pad = 0.1
        figure, tax, temp = density_dorito(temp, c,
                                 density_scale,
                                 density_cmap,
                                 density_vmax,
                                 log_density,
                                 pad=pad)
        scale = density_scale
        figure.set_size_inches(13,10)

    # if we're jittering, adjust the points for each thing
    if jitter:
        temp, c = jitter_dorito(temp, c, density_scale)

    # figure layout parameters
    fontsize = 18
    offset = 0.1
    mult = scale/5

    # if we don't already have a fig and axis from density,
    # make one
    if not density:
        figure, tax = ternary.figure(scale=scale, permutation='210')
        figure.set_facecolor('white')

    # plot gridlines below the scatterplot
    tax.gridlines(linewidth=3, multiple=mult,
                  color='white', zorder=1, linestyle=None)

    # scatter
    if scatter:
        figure, tax = scatter_dorito(temp, c, hue,
                                    size, log_size,
                                    cmap, mmap, alpha,
                                    density, legend,
                                    figure, tax)

    # sectors
    if sectors:
        line_dorito(sect_alpha, sect_beta, sect_gamma,
                    scale, tax, figure)

    # title handler
    if not title:
        if gene:
            title = '$\it{}$\n'.format(gene)
        else:
            title = ''
    else:
        if gene:
            title = '{} $\it{}$\n'.format(title, gene)
        else:
            title = '{}\n'.format(title)

    tax.set_title(title, fontsize=20)
    tax.boundary(linewidth=2, c='#e5ecf6')
    labels = ['{:.1f}'.format(n) for n in np.arange(0, 1.2, .2)]
    tax.ticks(ticks=labels,
              axis='lbr', linewidth=1, multiple=mult,
              tick_formats="%.1f", offset=0.014,
              fontsize=14)
    # tax.ticks(axis='lbr', linewidth=1, multiple=mult,
    #           tick_formats="%.1f", offset=0.014,
    #           fontsize=14)

    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    tax.set_background_color('#e5ecf6')

    if top == 'splicing_ratio':
        top_label = 'Splicing ratio $\\beta$'
    elif top == 'intron_chain':
        top_label = 'Intron chains $\\delta$'
    # tax.left_corner_label('# TSSs $\\alpha$', fontsize=fontsize)
    # tax.top_corner_label(top_label, fontsize=fontsize)
    # tax.right_corner_label('# TESs $\\gamma$', fontsize=fontsize)
    tax.left_axis_label('TSS $\\alpha$', fontsize=fontsize, offset=0.12)
    tax.right_axis_label(top_label, fontsize=fontsize, offset=0.12)
    tax.bottom_axis_label('TES $\\gamma$', fontsize=fontsize, offset=0.00)

    figure.set_facecolor('white')

    # tax.show()

    # save figure
    fname = opref
    if gene:
        fname += '_{}'.format(gene)
    if density:
        fname += '_density'
    if scatter:
        fname += '_scatter'
    if hue:
        fname += '_{}'.format(hue)
    fname += '.png'
    plt.savefig(fname, dpi=300, bbox_inches='tight')

    return temp
