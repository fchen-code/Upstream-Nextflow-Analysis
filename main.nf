#!/usr/bin/env nextflow
nextflow.enable.dsl=2
process fastqc {
	
	publishDir "fastqc_result", mode: "copy"

	input:
	path inputfile
	output:
	path "*fastqc.zip", emit: fastqc
	script:
	"""
	fastqc $inputfile
	"""
}

process multiqc {
	
	publishDir "multiqc_result", mode: "copy"
	
	input:
	path fastqc
	
	output:
	path "fastqc_result", type: "dir"
	
	script:
	"""
	mkdir fastqc_result
	multiqc . -o .
	"""
}

process trimming {
	publishDir "trimmed_result", mode: "copy"
	
	input:
	tuple val(name), path(file1), path(file2)
	
	output:
	path "*"
	path "*trimmed*", emit: "trimmed"

	script:
	"""
	fastp -i ${file1} -o ${name}_R1.trimmed.fastq -I ${file2} -O ${name}_R2.trimmed.fastq  
	"""
}

process star {
	
	publishDir "${projectDir}", mode: 'copy'
	
	input:
	tuple path(gtf), path(fa)
	
	output:
	path "*", type: 'dir', emit: index
	
	script:
	"""
	STAR --runThreadN 2 \\
	--runMode genomeGenerate \\
	--genomeDir STAR \\
	--genomeSAindexNbases 10 \\
	 --sjdbGTFfile ${gtf} \\
	--genomeFastaFiles ${fa}
	"""
}

process alignment {
	
	publishDir 'aligned_result', mode: 'copy'

	cpus 4
	
	input:
	tuple val(name), path(R1_file), path(R2_file), path(genomeDir)
	
	output:
	path '*'
	path '*.bam', emit:bams

	script:
	"""
	mkdir -p tmp
	chmod 700 tmp
	STAR --runMode alignReads \\
	--readFilesIn ${R1_file} ${R2_file} \\
	--genomeDir "${projectDir}/STAR" \\
	--outSAMtype BAM SortedByCoordinate \\
	--outFileNamePrefix "${name}" \\
	--readFilesCommand cat \\
	--outTmpDir "${projectDir}/tmp"
	"""
}

process bam_index {
	publishDir 'bam_index', mode: 'copy'
	
	input:
	path (file)

	output:
	path '*.bai', emit: bais

	script:
	"""
	samtools index ${file}
	"""
}

process feature_count {
	publishDir 'feature_counts', mode: 'copy'
	
	input:
	tuple path(gtf), path(bam)
	
	output:
	path '*'
	path '*txt.summary', emit: counts
	
	script:
	"""
	featureCounts -T 8 \\
	-p -B \\
	-s 0 \\
	-t exon \\
	-g gene_id \\
	-a ${gtf} \\
	-o count.txt \\
	${bam}
	"""
}

process final_qc {
	publishDir 'final_multiqc', mode: 'copy'
	
	input:
	tuple path(bam), path(count)
		
	output:
	path '*'
	
	script:
	"""
	multiqc ${bam} ${count}
	"""
}
params.data = "data/*.fastq"

params.gtf = 'gtf/*.gtf'
 
params.refs = 'refs/*.fa'

qc_ch = Channel.fromPath(params.data)

tr_ch2 = Channel.fromFilePairs('data/*_{1,2}.fastq', flat: true)
		.map{name, file1, file2 -> tuple(name, file1, file2)}

ge_ch = Channel.value(tuple(file(params.gtf), file(params.refs)))


workflow {

	qc = fastqc(qc_ch)
	multiqc(qc.fastqc.collect())
	trimmed = trimming(tr_ch2)
	ge_result = star(ge_ch)
	trimmed.trimmed.map
			{file1, file2 -> tuple("${file1.getFileName()}"
			.split("_R1")[0], file1, file2)}
			.combine(ge_result.index)
			.set{alignment_files}
	a_file = alignment(alignment_files)
	bam_index(a_file.bams)
	co_files = Channel.fromPath(params.gtf).map{file -> tuple(file)}.combine(a_file.bams)
	counts = feature_count(co_files)
	final_qc_files = a_file.bams.map{file -> tuple(file)}.combine(counts.counts)
	final_qc(final_qc_files)
}
