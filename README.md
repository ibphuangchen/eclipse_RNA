# This is a bulk RNAseq pipeline implantmented by snakemake and singularity

1. Create a working directory (WD) to be contain all the results (e.g. 'RNA_analysis')
2. Copy **Snakefile**, **cluster.jason**, and **snakemake.sh** into WD, and make sure the path of **RNAseq_config.ymal** is set correctly in **Snakefile**
3. Create a folder in WD to link the fastq raw files by "ln -s ...", and modify line 5 of **Snakefile** to ensure the right start fastq folder
4. run 
```bash
bio #for whatever reason the HPC numpy version is too old, so please enter the bio conda environment first 
```
```bash
sh snakemake.sh
```
 

