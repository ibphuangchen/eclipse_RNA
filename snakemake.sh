#!/usr/bin/env sh

set -e

[ ! -d bsub_log ] && mkdir -p bsub_log	

snakemake --use-singularity -j 500 --cluster-config ./cluster.jason \
					--rerun-incomplete \
					--singularity-args "-B /rsrch3,/home/p_eclipse_combio/workspace/,/home/p_eclipse_combio/.vep --no-home" \
		          --cluster 'bsub -n {cluster.n} -q {cluster.queue} -W {cluster.time} -M {cluster.memory} -R {cluster.resources} -o {cluster.output} -e {cluster.error}' $@

