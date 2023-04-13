import pandas as pd

df = pd.read_csv('config.tsv', '\t')

def get_df_col(wc, df, col):
    val = df.loc[df.dataset==wc.dataset, col].values[0]
    return val

files = df.fname.tolist()
samples = df.sample.tolist()

rule all:
    input:
        expand(config['proc']['fastq'],
               sample=samples)

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
