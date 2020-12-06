configfile: "config.json"

rule all:
	input: 
		expand("vcf_files/{genome}_{sample}.vcf.gz", genome=config["genomes"], sample=config["samples"])
		#		expand("bwa_aligned/{genome}_{sample}.sam", genome=glob_wildcards("ref_genome/{genome}.fasta").genome, sample=glob_wildcards("fastq_reads/{sample}_R1.fastq.gz").sample)



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

	threads: 4
	shell:
		"trimmomatic PE -threads 4 {input.fwd} {input.rev} {output.fwd_paired} {output.fwd_unpaired} {output.rev_paired} {output.rev_unpaired} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"

rule fastQC:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz"
	output:
		directory("fastqc/{sample}/")
	threads: 2
	conda:
		"envs/var_call.yaml"
	shell:
		"mkdir {output} && fastqc --quiet --threads 2 -o {output} -f fastq {input.fwd_paired} {input.rev_paired}"


rule bwa_align:
	input:
		fwd_paired="trimmomatic/{sample}_R1_paired.fastq.gz",
		rev_paired="trimmomatic/{sample}_R2_paired.fastq.gz",
		genome="ref_genome/{genome}.fasta",
		fastqc="fastqc/{sample}/"

	output:
		"bwa_aligned/{genome}_{sample}.sam"

	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"bwa mem -t 10 {input.genome} {input.fwd_paired} {input.rev_paired} > {output} && touch {input.fastqc}"

rule samtools_convert:
	input:
		"bwa_aligned/{genome}_{sample}.sam"

	output:
		"bwa_aligned/{genome}_{sample}.bam"

	threads: 1
	conda:
		"envs/var_call.yaml"

	shell:
		"samtools view -bS {input} > {output}"

rule samtools_fixmate:
	input:
		"bwa_aligned/{genome}_{sample}.bam"
	output:
		"bwa_aligned/{genome}_{sample}_fixmate.bam"

	threads: 5
	conda:
		"envs/var_call.yaml"
	shell:
		"samtools fixmate -@ 4 -O bam {output} {input}"

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

rule GATK_RealignerTargetCreator:
	input:
		genome="ref_genome/{genome}.fasta",
		bam="bwa_aligned/{genome}_{sample}_fixmate_sorted.bam"

	output:
		"gatk_realigner/{genome}_{sample}_realigner.intervals"
	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"GenomeAnalysisTK -T RealignerTargetCreator -R {input.genome} -I {input.bam} --num_threads 10 -o {output}"


rule GATK_IndelRealigner:
	input:
		genome="ref_genome/{genome}.fasta",
		bam="bwa_aligned/{genome}_{sample}_fixmate_sorted.bam",
		intervals="gatk_realigner/{genome}_{sample}_realigner.intervals"

	output:
		"gatk_realigner/{genome}_{sample}_realigned.bam"
	threads: 10
	conda:
		"envs/var_call.yaml"

	shell:
		"GenomeAnalysisTK -T IndelRealigner -R {input.genome} -I {input.bam} -targetIntervals {input.intervals} --num_treads 10 -o {output}"

rule samtools_indexing:
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
		ref_genome = "ref_genome/{genome}.fasta",
		bam = "gatk_realigner/{genome}_{sample}_realigned.bam"

	output:
		"vcf_files/{genome}_{sample}.vcf.gz"

	threads: 1
	conda:
		"envs/var_call.yaml"

	shell:
		"bcftools mpileup -Ou -f {input.ref_genome} {input.bam} | bcftools call -vmO z -o {output}"
