import os
import pandas as pd
logs='log'
configfile: "RNAseq_config.ymal"
fastqPath="bulkRNA/"
rnaSeqSamples = os.listdir(fastqPath)
rnaSeqSamples = pd.Series(rnaSeqSamples).str.replace('_L.*','').tolist()


def readStrand(infer_result):
	values = pd.read_table(infer_result, sep=':', skiprows=range(3), header=None)
	if values.iloc[1,1]/values.iloc[2,1]>5:
		return('FIRST_READ_TRANSCRIPTION_STRAND')
	elif values.iloc[1,1]<values.iloc[2,1]<0.2:
		return('SECOND_READ_TRANSCRIPTION_STRAND')
	else:
		return('NONE')

container: "/rsrch3/home/thera_dis/p_eclipse_combio/workspace/platform/containers/eclipse_rna_improve.sif"

rule all:
	input:
		"multiqc_out/multiqc_report.html"
		
rule combineFastq:
	input:
		R1L1=fastqPath+"{sample}_L001_R1_001.fastq.gz",
		R1L2=fastqPath+"{sample}_L002_R1_001.fastq.gz",
		R2L1=fastqPath+"{sample}_L001_R2_001.fastq.gz",
		R2L2=fastqPath+"{sample}_L002_R2_001.fastq.gz"
	output:
		R1=temp("fastq/{sample}.R1.fastq.gz"),
		R2=temp("fastq/{sample}.R2.fastq.gz")
	shell:
		"""
		ls -l /home &> debug.log
		ls -l /home/p_eclipse_combio &>> debug.log
		ls -l /home/p_eclipse_combio/workspace/runs/220429_A00667_0052_AHC75LDRX2/fastq/bulkRNA/H2595_S15_L001_R1_001.fastq.gz &>> debug.log
		cat {input.R1L1} {input.R1L2} > {output.R1}
		cat {input.R2L1} {input.R2L2} > {output.R2}
		"""

rule fastqc:
	input:
		"fastq/{sample}.R1.fastq.gz",
		"fastq/{sample}.R2.fastq.gz",
		"processed_fastq/{sample}_R1.fastq.gz",
		"processed_fastq/{sample}_R2.fastq.gz"
	threads:
		12
	output:
		"fastqc/{sample}.R1_fastqc.html",
		"fastqc/{sample}_R1_fastqc.html"
	log: "log/{sample}_fastqc.log" 
	shell:
		"fastqc -t {threads} {input} -o fastqc 2>&1>{log}"

rule bbduk_process:
	input:
		fastq1="fastq/{sample}.R1.fastq.gz", 
		fastq2="fastq/{sample}.R2.fastq.gz"
	output:
		first=temp("processed_fastq/{sample}_R1.fastq.gz"),
		second=temp("processed_fastq/{sample}_R2.fastq.gz"),
		stats='bbmap/{sample}_stats.txt'
	log: "log/{sample}_bbduk.log"
	threads:
		12
	shell:
		"bbduk.sh -Xmx20g -Xms4g t=12 overwrite=t \
        in1={input.fastq1} in2={input.fastq2} \
        out1={output.first} out2={output.second} \
		ref=/rsrch3/home/thera_dis/p_eclipse_combio/.linuxbrew/opt/bbtools/resources/adapters.fa \
        ktrim=r k=23 mink=11 hdist=1 entropy=0.2 entropywindow=50 stats={output.stats} tbo 2>&1>{log}"

rule optitype:
	input:
		fq1="processed_fastq/{sample}_R1.fastq.gz",
		fq2="processed_fastq/{sample}_R2.fastq.gz"
	output:
		"hlatype/{sample}_result.tsv"
	singularity:
		"/home/p_eclipse_combio/workspace/platform/containers/optitype_latest.sif"
	shell:
		"""
		OptiTypePipeline.py -i {input.fq1} {input.fq2} \
		-o hlatype -r --prefix {wildcards.sample}
		"""

rule star_align:
	input:
		fastq1='processed_fastq/{sample}_R1.fastq.gz',
		fastq2='processed_fastq/{sample}_R2.fastq.gz'
	output:
		temp("starAlign/{sample}_Aligned.sortedByCoord.out.bam"),
		temp("starAlign/{sample}_Aligned.toTranscriptome.out.bam")
	threads:
		12
	log: "log/{sample}_starAlign.log"
	shell:
		"STAR \
		--runMode alignReads \
        --runThreadN {threads} \
        --genomeDir {config[starRef]} \
        --readFilesIn {input.fastq1} {input.fastq2}\
        --readFilesCommand zcat \
        --outSAMunmapped Within KeepPairs \
        --twopassMode Basic \
		--outSAMtype BAM SortedByCoordinate\
		--quantMode TranscriptomeSAM \
		--chimOutType WithinBAM SoftClip --chimJunctionOverhangMin 20 --chimScoreMin 1 --chimScoreDropMax 30 \
		--chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --alignSJstitchMismatchNmax 5 -1 5 5 --chimSegmentReadGapMax 3 \
		--chimSegmentMin 10 \
        --outFileNamePrefix starAlign/{wildcards.sample}_ 2>&1>{log}"

rule samIndex:
	input:
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam" 
	output:
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam.bai"
	shell:
		"samtools index {input}"

rule arriba:
	input:
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		keep='fusion/{sample}_arriba_fusions.tsv',
		discard='fusion/{sample}_arriba_fusions.discarded.tsv'
	log: "log/{sample}_arriba.log"
	shell:
		"arriba -x {input} -o {output.keep} -O {output.discard} \
		-a /rsrch3/scratch/thera_dis/chuang8/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		-g {config[ensemblGTF]} \
		-b {config[arribaBlack]} 2>&1>{log}"

rule filterUmmapedSAM:
	input:
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		temp("starAlign/{sample}_Aligned.sortedByCoord.out.filter.bam")
	threads:
		8
	shell:
		"samtools view -@ {threads} -b -F 4 {input} > {output}"


rule inferStrand:
	input:
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam",
		"starAlign/{sample}_Aligned.sortedByCoord.out.bam.bai"
	output:
		"starAlign/{sample}_infer_experiment.log"
	shell:
		"infer_experiment.py -i {input} -s 10000000 -r {config[ensemblBED]} > {output}"

rule picardCollectRnaMetric:
	input:
		bam="starAlign/{sample}_Aligned.sortedByCoord.out.filter.bam",
		strand='starAlign/{sample}_infer_experiment.log'
	output:
		"picardCollectRnaMeric/{sample}.RNA_Metrics"
	log: "log/{sample}_picardCollectRNAMetric.log"
	run:
		strandParam=readStrand(input['strand'])
		shell("java -Xmx4g -jar {config[picard]} CollectRnaSeqMetrics \
		I={input.bam} O={output} \
		REF_FLAT={config[picardref]} \
		STRAND={strandParam} \
		RIBOSOMAL_INTERVALS={config[picardribo]}")

rule plotFusion:
	input:
		fusion="fusion/{sample}_arriba_fusions.tsv",
		bam="starAlign/{sample}_Aligned.sortedByCoord.out.bam"
	output:
		"fusion/{sample}_arriba_fusions.pdf"
	shell:
		"draw_fusions.R --fusions={input.fusion} \
		--alignments={input.bam} \
		--output={output} \
		--annotation={config[ensemblGTF]} \
		--cytobands={config[cytobands]} \
		--proteinDomains={config[proteinDomains]}"

rule RSEM:
	input:
		"starAlign/{sample}_Aligned.toTranscriptome.out.bam"
	output:
		"RSEM/{sample}.genes.results"
	threads: 8
	shell:
		"rsem-calculate-expression -p {threads} \
						--estimate-rspd \
						--append-names \
						--no-bam-output \
						--alignments --paired-end {input} \
						{config[rsemRef]} \
						RSEM/{wildcards.sample}"


rule callSnpsGenome:
	input:
		bam="starAlign/{sample}_Aligned.sortedByCoord.out.bam",
		ref_fa=config["humanGenomeFa"]
	output:
		snp=temp("snp/{sample}.snp.genome.vcf"),
		indel=temp("snp/{sample}.indel.genome.vcf")
	message: 
		"Running varscan for snp analysis genome wide"
	benchmark:
		"benchmarks/{sample}.call_snps_genome.txt"
	shell:
		"""samtools mpileup -f {input.ref_fa} {input.bam} | awk \'$4 != 0\' | \
		tee >(varscan mpileup2snp {wildcards.sample}.pileup --min-coverage 20 --min-reads2 4 --output-vcf > {output.snp}) | \
		varscan mpileup2indel {wildcards.sample}.pileup --min-coverage 20 --min-reads2 4 --output-vcf > {output.indel}
		"""

rule mergeVCF:
	input:
		snp="snp/{sample}.snp.genome.vcf",
		indel="snp/{sample}.indel.genome.vcf"
	output:
		"snp/{sample}.genome.vcf"
	shell:
		"""
		bgzip {input.snp}
		bgzip {input.indel}
		tabix -p vcf {input.snp}.gz
		tabix -p vcf {input.indel}.gz
		vcf-merge {input.snp}.gz {input.indel}.gz > {output}
		"""

rule vep:
	input:
		"snp/{sample}.genome.vcf"
	output:
		"snp/{sample}.genome.vep.vcf"
	singularity:
		"/rsrch3/scratch/thera_dis/chuang8/shared/singularity-images/ensembl-vep.simg"
	shell:
		"""
		vep --species homo_sapiens --assembly GRCh38 \
			--no_progress --no_stats --buffer_size 5000 \
			--sift b --ccds --uniprot --hgvs --symbol \
			--numbers --domains --gene_phenotype --canonical \
					--protein --biotype --uniprot --tsl \
					--variant_class --shift_hgvs 1 --check_existing \
					--total_length --allele_number --no_escape --xref_refseq \
					--failed 1 --vcf --flag_pick_allele \
					--pick_order canonical,tsl,biotype,rank,ccds,length \
					--dir /home/p_eclipse_combio/.vep \
					--fasta /home/p_eclipse_combio/.vep/homo_sapiens/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
					--format vcf --input_file {input} --output_file {output} \
					--offline --pubmed --fork 4 --cache_version 104 \
					--polyphen b --af --af_1kg --af_esp --af_gnomad --regulatory
			"""

rule vcf2vep:
	input:
		"snp/{sample}.genome.vep.vcf"
	output:
		"maf/{sample}.maf.tsv"
	shell:
		"""
		vcf2maf.pl --input-vcf {input} \
		--output-maf {output} \
		--inhibit-vep \
		--ncbi-build {config[vepRefBuild]} \
		--vep-data {config[vepCache]} \
		--species homo_sapiens \
		--ref-fasta {config[vepRefFasta]} \
		--tumor-id Sample1 \
		"""

rule multiqc:
	input:
		expand('bbmap/{s}_stats.txt', s=rnaSeqSamples),
		expand("fusion/{s}_arriba_fusions.tsv",s=rnaSeqSamples),
		expand("fusion/{s}_arriba_fusions.pdf",s=rnaSeqSamples),
		expand("fastqc/{s}.R1_fastqc.html",s=rnaSeqSamples),
		expand("fastqc/{s}_R1_fastqc.html",s=rnaSeqSamples),
		expand("picardCollectRnaMeric/{s}.RNA_Metrics",s=rnaSeqSamples),
		expand("starAlign/{s}_infer_experiment.log", s=rnaSeqSamples),
		expand("RSEM/{s}.genes.results", s=rnaSeqSamples),
		expand("hlatype/{s}_result.tsv", s=rnaSeqSamples),
		expand("maf/{s}.maf.tsv", s=rnaSeqSamples)	
	output:
		"multiqc_out/{sample}.html"
	threads:
		1
	resources:
		mem_mb=1000
	shell:
		"""
		multiqc . -f -o multiqc_out
		"""

