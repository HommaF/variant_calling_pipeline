configfile: "config.json"

rule all:
	input: 
		expand("vcf_files/{genome}_{sample}_filtered_multiallelic.bcf", genome=config["genomes"], sample=config["samples"])



rule trimmomatic:
	input:
		fwd="fastq_reads/{sample}_R1.fastq.gz",
		rev="fastq_reads/{sample}_R2.fastq.gz"
	output:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		fwd_unpaired="trimmomatic/{sample}_R1_unpaired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz",
		rev_unpaired="trimmomatic/{sample}_R2_unpaired.fastq.gz"

	conda:
		"envs/snp_calling.yaml"

	threads: 12
	shell:
		"trimmomatic PE -threads 12 {input.fwd} {input.rev} {output.fwd_paired} {output.fwd_unpaired} {output.rev_paired} {output.rev_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule fastQC:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		directory("fastqc/{sample}")
	threads: 2
	conda:
		"envs/var_call.yaml"
	shell:
		"mkdir {output} && fastqc --quiet --threads 2 -o {output} -f fastq {input.fwd_paired} {input.rev_paired}"


rule bwa_align:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz",
		genome="ref_genome/{genome}.fa",
		fastqc="fastqc/{sample}/"

	output:
		"bwa_aligned/{genome}_{sample}.sam"

	params:
		rg="@RG\\tID:{sample}\\tSM:{sample}"


	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"bwa mem -R '{params.rg}' -t 10 {input.genome} {input.fwd_paired} {input.rev_paired} > {output}"


rule samtools_fixmate:
	input:
		"bwa_aligned/{genome}_{sample}.sam"
	output:
		"bwa_aligned/{genome}_{sample}_fixmate.bam"

	threads: 10
	conda:
		"envs/var_call.yaml"
	shell:
		"samtools fixmate -@ 9 -O bam {input} {output}"

rule sam_compress:
	input:
		bam="bwa_aligned/{genome}_{sample}_fixmate.bam",
		sam="bwa_aligned/{genome}_{sample}.sam"

	output:
		"bwa_aligned/{genome}_{sample}.sam.gz"

	threads:1
	conda:
		"envs/var_call.yaml"
	shell:
		"bgzip {input.sam}"


rule samtools_sorting:
	input:
		"bwa_aligned/{genome}_{sample}_fixmate.bam"
	output:
		"bwa_aligned/{genome}_{sample}_fixmate_sorted.bam"

	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"samtools sort -@ 9 -O bam -o {output} {input}"


rule samtools_indexing_one:
	input:
		"bwa_aligned/{genome}_{sample}_fixmate_sorted.bam"
	output:
		"bwa_aligned/{genome}_{sample}_fixmate_sorted.bam.bai"

	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"samtools index -@ 10 {input}"


rule GATK_RealignerTargetCreator:
	input:
		genome="ref_genome/{genome}.fa",
		bam="bwa_aligned/{genome}_{sample}_fixmate_sorted.bam",
		bam_sorted="bwa_aligned/{genome}_{sample}_fixmate_sorted.bam.bai"

	output:
		"gatk_realigner/{genome}_{sample}_realigner.intervals"
	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"GenomeAnalysisTK -T RealignerTargetCreator -R {input.genome} -I {input.bam} --num_threads 10 -o {output}"


rule GATK_IndelRealigner:
	input:
		genome="ref_genome/{genome}.fa",
		bam="bwa_aligned/{genome}_{sample}_fixmate_sorted.bam",
		intervals="gatk_realigner/{genome}_{sample}_realigner.intervals"

	output:
		"gatk_realigner/{genome}_{sample}_realigned.bam"
	threads: 1
	conda:
		"envs/var_call.yaml"

	shell:
		"GenomeAnalysisTK -T IndelRealigner -R {input.genome} -I {input.bam} -targetIntervals {input.intervals} -o {output}"


rule samtools_indexing_two:
	input:
		"gatk_realigner/{genome}_{sample}_realigned.bam"
	output:
		"gatk_realigner/{genome}_{sample}_realigned.bam.bai"
	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"samtools index -@ 9 {input}"

rule bcf_mpileup:
	input:
		ref_genome = "ref_genome/{genome}.fa",
		bam = "gatk_realigner/{genome}_{sample}_realigned.bam",
		sam_comp="bwa_aligned/{genome}_{sample}.sam.gz"

	output:
		"vcf_files/{genome}_{sample}.bcf"

	threads: 2
	conda:
		"envs/var_call.yaml"

	shell:
		"bcftools mpileup --threads 2 -Ob -f {input.ref_genome} {input.bam} -o {output}"	

rule bcf_call:
	input:
		"vcf_files/{genome}_{sample}.bcf"
	output:
		"vcf_files/{genome}_{sample}_filtered.bcf"
	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"bcftools call --threads 10 -vcOb {input} | bcftools filter --threads 10 -Ob -o {output} -s LOWQUAL -i'%QUAL>20'"

rule bcf_norm:
	input:
		"vcf_files/{genome}_{sample}_filtered.bcf"
	output:
		"vcf_files/{genome}_{sample}_filtered_multiallelic.bcf"
	threads: 1
	conda:
		"envs/var_call.yaml"

	shell:
		"bcftools norm -m +any -Ob {input} > {output}"
