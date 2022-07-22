# A bulk RNAseq pipeline implantmented by snakemake and singularity

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

## Ouput
1. fastqc
2. RNA fustion (by Arriba)
3. Gene expression (by RSEM)
4. HLA typing (by Optitype)
5. Mutation (germine and somatic mixed, called by VarScan2 and annotated by vep and  vcf2maf)
6. Picard QC
7. MultiQC report

 

