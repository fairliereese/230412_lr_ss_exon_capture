import pandas as pd

configfile: 'config.yml'

df = pd.read_csv('config.tsv', '\t')

def get_df_col(wc, df, col):
    val = df.loc[df.dataset==wc.dataset, col].values[0]
    return val

files = df.fname.tolist()
samples = df['sample'].tolist()
datasets = df['dataset'].tolist()
platforms = df['platform'].tolist()

rule all:
    input:
        # expand(config['proc']['sam_tag'],
        #        dataset=datasets),
        config['proc']['demux_db'],
        # expand(config['proc']['sam'],
        #       dataset=datasets),
        config['proc']['db']

rule symlink:
    resources:
        mem_gb = 4,
        threads = 1
    shell:
        "ln -s {params.fname} {output.out}"

use rule symlink as sl_fastq with:
    params:
      fname = lambda wc:get_df_col(wc, df, 'fname')
    output:
      out = config['proc']['fastq']

rule demux:
  input:
      fq = config['proc']['fastq']
  output:
      dx_fq = config['proc']['demux_fastq']
  resources:
      mem_gb = 32,
      threads = 4
  params:
      d = config['lr_splitpipe'],
      opref = config['proc']['demux_fastq'].rsplit('_demux.fastq', maxsplit=1)[0]
  shell: """python {params.d}demultiplex.py all \
      -f {input.fq} \
      -o {params.opref} \
      -t {resources.threads} \
      -k WT_mega \
      --l1_mm 4 \
      --l2_mm 4 \
      --max_read_len 10000 \
      --max_linker_dist 200 \
      --chunksize 2000000 \
      --verbosity 2 \
      --delete_input
  """

  # mapping - versions w/ and w/o demux just to get things
  # started a lil earlier

rule map:
  resources:
    threads = 32,
    mem_gb = 64
  shell:
      'minimap2 --MD \
     			-t {resources.threads} \
     			-ax splice \
     			-k14 \
     		    {input.ref_fa} {input.fastq} > {output.sam}'

use rule map as map_demux with:
    input:
        fastq = config['proc']['demux_fastq'],
        ref_fa = config['ref']['fa']
    output:
        sam = config['proc']['demux_sam']

use rule map as map_not_demux with:
    input:
        fastq = config['proc']['fastq'],
        ref_fa = config['ref']['fa']
    output:
        sam = config['proc']['sam']

# replace bam tag
rule bam_tag:
    input:
        sam = config['proc']['demux_sam']
    resources:
        threads = 8,
        mem_gb = 32
    params:
        d = config['lr_splitpipe'],
        opref = config['proc']['sam_tag'].replace('.sam', '')
    output:
        sam = config['proc']['sam_tag']
    shell:
        """python {params.d}add_bam_tag.py \
            -s {input.sam} \
            -k WT_mega \
            --merge_primers \
            --suffix {wildcards.dataset} \
            -o
        """

# talon - versions w/ and w/o demux just to get things
# started a lil earlier
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
        sam_files = expand(config['proc']['sam_tag'], dataset=datasets)
    resources:
        threads = 1,
        mem_gb = 1
    params:
        samples = samples,
        platforms = platforms
    output:
        talon_config = config['proc']['demux_config']
    run:
        dumb = [[s,p,f] for s,p,f in zip(params.samples,
                                         params.platforms,
                                         params.sam_files)]
        cfg_df = pd.DataFrame(data=dumb)
        cfg_df.to_csv(talon_config, sep=',', index=False, header=None)

rule talon_config:
    input:
        sam_files = expand(config['proc']['sam'], dataset=datasets)
    params:
        datasets = datasets,
        samples = samples,
        platforms = platforms
    resources:
        threads = 1,
        mem_gb = 1
    output:
        talon_config = config['proc']['config']
    run:
        dumb = [[d,s,p,f] for d,s,p,f in zip(params.dataset,
                                         params.samples,
                                         params.platforms,
                                         params.sam_files)]
        cfg_df = pd.DataFrame(data=dumb)
        cfg_df.to_csv(talon_config, sep=',', index=False, header=None)


rule talon_cb:
    resources:
        threads = 32,
        mem_gb = 256
    shell:
        "talon \
          --f {input.cfg} \
          --cb \
          --db {input.annot_db} \
          --build {params.build} \
          -t {resources.threads} \
          -c 0.8 \
          --o {output.db}"

use rule talon_cb as talon_demux with:
    input:
        annot_db = config['proc']['init_db'],
        cfg = config['proc']['demux_config']
    params:
        build = 'mm10',
        opref = config['proc']['demux_db'].replace('.db', '')
    output:
        db = config['proc']['demux_db']

rule talon:
    resources:
        threads = 32,
        mem_gb = 256
    shell:
        "talon \
          --f {input.cfg} \
          --db {input.annot_db} \
          --build {params.build} \
          -t {resources.threads} \
          -c 0.8 \
          --o {output.db}"

use rule talon as talon_not_demux with:
  input:
      annot_db = config['proc']['init_db'],
      cfg = config['proc']['config']
  params:
      build = 'mm10',
      opref = config['proc']['db'].replace('.db', '')
  output:
      db = config['proc']['db']
