```bash
snakemake \
  -s Snakefile \
  -j 20 \
  --latency-wait 120 \
  --cluster-cancel scancel \
  --cluster "sbatch -A seyedam_lab --partition=standard --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n

snakemake \
  -s Snakefile \
  -j 20 \
  --latency-wait 120 \
  --cluster "sbatch -A seyedam_lab --partition=highmem --mem={resources.mem_gb}GB -c {resources.threads} --mail-user=freese@uci.edu --mail-type=START,END,FAIL --time=72:00:00" -n


snakemake \
  -s Snakefile \
  -j 10 \
  --latency-wait 120 \
  -n
```
