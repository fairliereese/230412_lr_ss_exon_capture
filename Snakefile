import pandas as pd
import os
import sys
import re

# p = os.path.dirname(os.getcwd())
p = os.getcwd()
sys.path.append(p)

from utils import *

configfile: 'config.yml'
config_tsv = '230427_config.tsv'
auto_dedupe = True

df = parse_config_file(config_tsv, auto_dedupe=auto_dedupe)

# subset the config df on dataset
def get_df_dataset(dataset, df):
    return df.loc[df.dataset==dataset]

# subset the config df on flowcell
def get_df_flowcell(flowcell, df):
    return df.loc[df.flowcell==flowcell]

# get a 1:1 value for dataset:<col> from config
def get_df_dataset_val(wc, df, col):
    try:
        dataset = wc.dataset
    except:
        dataset = wc
    temp = get_df_dataset(dataset, df)
    return temp[col].values[0]

# get a 1:many value (ie dataset:flowcell, dataset:fname) from config
def get_df_dataset_col(wc, df, col):
    try:
        dataset = wc.dataset
    except:
        dataset = wc
    temp = get_df_dataset(dataset, df)
    return temp[col].tolist()

# get a 1:1 value from flowcell,dataset:<col>
def get_df_dataset_flowcell_col(wc, df, col):
    temp = get_df_dataset(wc.dataset, df)
    temp = get_df_flowcell(wc.flowcell, temp)
    return temp[col].values[0]

# def get_df_col(wc, df, col):
#     val = df.loc[df.dataset==wc.dataset, col].values[0]
#     return val
#
# def get_df_whole_col(wc, df, col):
#     temp = df.loc[df.dataset==wc.dataset]
#     if flowcell in wc.keys():
#         temp = df.loc[df.flowcell==wc.flowcell]
#     vals = temp[col].tolist()
#     return vals

def get_sublib_bc_files(wc, df, config, str=False):
    sublib_flowcells = get_df_dataset_col(wc, df, 'flowcell')
    bc_files = expand(config['proc']['demux_bc'],
                      dataset=wc.dataset,
                      flowcell=sublib_flowcells)
    if str:
        bc_files = ','.join(bc_files)
    return bc_files

files = df['fname'].tolist()
samples = df['sample'].unique().tolist()
datasets = df['dataset'].unique().tolist()
platforms = df['platform'].unique().tolist()
flowcells = df['flowcell'].tolist()

wildcard_constraints:
    dataset= '|'.join([re.escape(x) for x in datasets])

rule all:
    input:
        expand(config['proc']['map_stats'],
              dataset=datasets),
        expand(config['proc']['sam_rev_stats'],
              dataset=datasets),
        expand(config['proc']['tc_stats'],
              dataset=datasets),
        config['proc']['adata'],
        config['proc']['g_adata']

# rule symlink:
#     resources:
#         mem_gb = 4,
#         threads = 1
#     shell:
#         "ln -s {params.fname} {output.out}"
#
# use rule symlink as sl_fastq with:
#     params:
#       fname = lambda wc:get_df_dataset_flowcell_col(wc, df, 'fname')
#     output:
#       out = config['proc']['fastq']

################################################################################
############################# LR-Splitpipe #####################################
################################################################################


# rule demux:
#   # input:
#   #     fq = config['proc']['fastq']
#   output:
#       dx_fq = config['proc']['demux_fastq']
#   resources:
#       mem_gb = 32,
#       threads = 4
#   params:
#       fname = lambda wc:get_df_dataset_flowcell_col(wc, df, 'fname'),
#       d = config['lr_splitpipe'],
#       opref = config['proc']['demux_fastq'].rsplit('_demux.fastq', maxsplit=1)[0]
#   shell: """python {params.d}demultiplex.py all \
#       -f {input.fq} \
#       -o {params.opref} \
#       -t {resources.threads} \
#       -k WT \
#       -c v2 \
#       --l1_mm 4 \
#       --l2_mm 4 \
#       --max_read_len 10000 \
#       --max_linker_dist 200 \
#       --chunksize 2000000 \
#       --verbosity 2 \
#       --delete_input
#   """

rule demux_find_bcs:
  # input:
  #     fq = config['proc']['fastq']
  output:
      bc = config['proc']['demux_bc']
  resources:
      mem_gb = 128,
      threads = 16
  params:
      fname = lambda wc:get_df_dataset_flowcell_col(wc, df, 'fname'),
      d = config['lr_splitpipe'],
      opref = config['proc']['demux_bc'].rsplit('_bcs.tsv', maxsplit=1)[0]
  shell:
      """python {params.d}demultiplex.py find_bcs \
          -f {params.fname} \
          -o {params.opref} \
          -t {resources.threads} \
          -k WT \
          -c v2 \
          --l1_mm 4 \
          --l2_mm 4 \
          --max_read_len 10000 \
          --max_linker_dist 200 \
          --chunksize 2000000 \
          --verbosity 2 \
          --delete_input
      """

rule demux_proc_bcs:
    input:
        bc_files = lambda wc: get_sublib_bc_files(wc, df, config)
    output:
        fq = config['proc']['demux_fastq']
    resources:
        mem_gb = 64,
        threads = 4
    params:
        d = config['lr_splitpipe'],
        opref = config['proc']['demux_fastq'].rsplit('_demux.fastq', maxsplit=1)[0],
        str_bc_files = lambda wc: get_sublib_bc_files(wc, df, config, str=True)
    shell:
        """
        python {params.d}demultiplex.py process_bcs \
            -f {params.str_bc_files} \
            -o {params.opref} \
            -t {resources.threads} \
            -k WT \
            -c v2 \
            --chunksize 2000000 \
            --verbosity 2 \
            --delete_input
        """

################################################################################
################################ Mapping #######################################
################################################################################

# mapping - versions w/ and w/o demux just to get things
# started a lil earlier
rule map:
  resources:
    threads = 16,
    mem_gb = 32
  shell:
      """
      module load minimap2
      minimap2 --MD \
               -t {resources.threads} \
               -ax splice \
               -k14 \
               {input.ref_fa} {input.fastq} > {output.sam} 2> {output.log}"""

rule alignment_stats:
   resources:
       threads = 2,
       mem_gb = 16
   shell:
       """
       module load samtools
       samtools stats {input.alignment} | grep ^SN | cut -f 2- | grep -e 'reads map
ped' -e 'reads unmapped' -e 'average length' -e 'maximum length' | sed '/reads mapped and paired/d' > {output.stats}
       """

rule sam_to_bam:
    resources:
        threads = 16,
        mem_gb = 16
    shell:
        """
        module load samtools
        samtools sort \
            --threads {params.threads} \
            -O bam {input.sam} > {output.bam}
        samtools index -$ {params.threads} {output.bam}
        """

use rule map as map_demux with:
    input:
        fastq = config['proc']['demux_fastq'],
        ref_fa = config['ref']['fa']
    output:
        sam = config['proc']['sam'],
        log = config['proc']['sam_log']

use rule alignment_stats as map_stats with:
    input:
        alignment = config['proc']['sam']
    output:
        stats = config['proc']['map_stats']

rule rev_alignment:
    input:
        sam = config['proc']['sam']
    resources:
        threads = 8,
        mem_gb = 32
    output:
        sam_rev = config['proc']['sam_rev']
    run:
        reverse_alignment(input.sam, output.sam_rev, resources.threads)

use rule alignment_stats as rev_sam_stats with:
    input:
        alignment = config['proc']['sam_rev']
    output:
        stats = config['proc']['sam_rev_stats']


################################################################################
############################# TranscriptClean ##################################
################################################################################

rule tc:
    resources:
        mem_gb = 80,
        threads = 16
    shell:
        """
        python {params.tc}TranscriptClean.py \
            -t {resources.threads} \
            --sam {input.sam} \
            --genome {input.fa} \
            --correctIndels True \
            --canonOnly \
            --primaryOnly \
            --deleteTmp \
            --tmpDir {params.opref}_temp/ \
            --outprefix {params.opref} 2> {output.log}
        """

use rule tc as tc_sam with:
    input:
        sam = config['proc']['sam_rev'],
        fa = config['ref']['fa']
    params:
        tc = config['tc_path'],
        opref = config['proc']['sam_clean'].rsplit('_clean.sam', maxsplit=1)[0]
    output:
        sam = config['proc']['sam_clean'],
        log = config['proc']['tc_log']

use rule alignment_stats as tc_stats with:
    input:
        alignment = config['proc']['sam_clean']
    output:
        stats = config['proc']['tc_stats']

################################################################################
############################### CB tag BAM #####################################
################################################################################

# replace bam tag
rule bam_tag:
    input:
        sam = config['proc']['sam_clean']
    resources:
        threads = 8,
        mem_gb = 32
    params:
        d = config['lr_splitpipe'],
        opref = config['proc']['sam_tag'].replace('_merged_primers.sam', '')
    output:
        sam = config['proc']['sam_tag']
    shell:
        """python {params.d}add_bam_tag.py \
            -s {input.sam} \
            -k WT \
            -c v2 \
            --merge_primers \
            --suffix {wildcards.dataset} \
            -o {params.opref}
        """

################################################################################
################################# TALON ########################################
################################################################################
rule talon_label:
    input:
        fa = config['ref']['fa'],
        sam = config['proc']['sam_tag']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        opref = config['proc']['sam_label'].rsplit('_labeled.sam', maxsplit=1)[0]
    output:
        sam = config['proc']['sam_label']
    shell:
        """
        talon_label_reads \
            --f {input.sam} \
            --g {input.fa} \
            --tmpDir {params.opref} \
            --ar 20 \
            --deleteTmp \
            --o {params.opref}
        """

rule talon_init:
    input:
        ref_gtf = config['ref']['gtf']
    resources:
        threads = 16,
        mem_gb = 64
    params:
        build = 'mm10',
        annot_ver = 'vM21',
        opref = config['proc']['init_db'].replace('.db', '')
    output:
        out = config['proc']['init_db']
    shell:
        "talon_initialize_database \
        --f {input.ref_gtf} \
        --g {params.build} \
        --a {params.annot_ver} \
        --l 0 \
        --idprefix lr_splitseq \
        --5p 500 \
        --3p 300 \
        --o {params.opref}"

rule talon_cb_config:
    input:
        sam_files = expand(config['proc']['sam_label'], dataset=datasets)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        # samples = samples,
        # platforms = platforms
        # df = df
        # samples = lambda wc:get_df_dataset_val(wc, df, 'sample'),
        # platforms = lambda wc:get_df_dataset_val(wc, df, 'platform')
        datasets = datasets
    output:
        talon_config = config['proc']['demux_config']
    run:
        # dumb = [[s,p,f] for s,p,f in zip(params.samples,
        #                                  params.platforms,
        #                                  input.sam_files)]
        # cfg_df = pd.DataFrame(data=dumb)
        samples = []
        platforms = []
        for d in datasets:
            samples.append(get_df_dataset_val(d, df, 'sample'))
            platforms.append(get_df_dataset_val(d, df, 'platform'))
        zipped = [[s,p,f] for s,p,f in zip(samples,
                                           platforms,
                                           input.sam_files)]
        cfg_df = pd.DataFrame(data=zipped)
        cfg_df.to_csv(output.talon_config, sep=',', index=False, header=None)

rule talon_cb:
    resources:
        threads = 16,
        mem_gb = 360
    shell:
        """
        # copy the input db so that a copy of it will remain unmodified
        cp {input.annot_db} {input.annot_db}_back

        talon \
          --f {input.cfg} \
          --cb \
          --db {input.annot_db} \
          --build {params.build} \
          -t {resources.threads} \
          -c 0.8 \
          --o {params.opref} \
          --tmpDir {params.opref}_tmp

          # mv the output db with annotations
          # and restore the annot db
          mv {input.annot_db} {output.db}
          mv {input.annot_db}_back {input.annot_db}"""

use rule talon_cb as talon_demux with:
    input:
        annot_db = config['proc']['init_db'],
        cfg = config['proc']['demux_config']
    params:
        build = 'mm10',
        opref = config['proc']['demux_db'].rsplit('.db', maxsplit=1)[0]
    output:
        db = config['proc']['demux_db']

# rule talon_filt:
#   input:
#       db = config['proc']['demux_db']
#   resources:
#       threads = 1,
#       mem_gb = 32
#   params:
#       annot = 'vM21'
#   output:
#       list = config['proc']['filt_list']
#   shell:
#       """
#       talon_filter_transcripts \
#           --db {input.db} \
#           -a {params.annot} \
#           --maxFracA=0.5 \
#           --minCount=1 \
#           --minDatasets=4 \
#           --o {output.list}
#       """

rule talon_adata:
    input:
        db = config['proc']['demux_db'],
        # pass_list = config['proc']['filt_list']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        annot = 'vM21',
        build = 'mm10'
    output:
        h5ad = config['proc']['adata']
    shell:
        """
        talon_create_adata \
            --db {input.db} \
            -a {params.annot} \
            -b {params.build} \
            --o {output.h5ad}
        """

rule talon_gene_adata:
    input:
        db = config['proc']['demux_db'],
        # pass_list = config['proc']['filt_list']
    resources:
        threads = 1,
        mem_gb = 32
    params:
        annot = 'vM21',
        build = 'mm10'
    output:
        h5ad = config['proc']['g_adata']
    shell:
        """
        talon_create_adata \
            --db {input.db} \
            -a {params.annot} \
            -b {params.build} \
            --gene \
            --o {output.h5ad}
        """
