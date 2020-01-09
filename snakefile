from glob import glob

def seq_names():
    return [ '.'.join(path.split('/')[-1].split('.')[:-1])
		for path in glob(checkpoints.fetch_source.get().output[0] + '/*.fastq')
		]

def bams(_):
    seqs = seq_names()
    return [
		'mapped/%s.Aligned.sortedByCoord.out.bam' % seq[:-3] for seq
                in seqs
		]

rule collect:
	input:
		'counts.txt'

checkpoint fetch_source:
	output:
		directory('source')
	shell:
		"""
		mkdir -p {output}
		curl 'http://rssf.i2bc.paris-saclay.fr/X-fer/AtelierNGS/TPrnaseq.tar.gz' | tar xzf - -C {output}
		"""

rule fastqc:
	input:
		'{seq}.fastq'
	output:
		'quality/{seq}_fastqc.html'
	shell:
		"""
		fastqc "{input}" -o `dirname {output}`
		"""

rule trim:
	input:
		'source/{seq}.R1.fastq',
		'source/{seq}.R2.fastq'
	output:
		'trimmed/{seq}.R1P.fastq',
		'trimmed/{seq}.R1U.fastq',
		'trimmed/{seq}.R2P.fastq',
		'trimmed/{seq}.R2U.fastq'
	shell:
		"""
		trimmomatic PE {input} {output} LEADING:20 TRAILING:20 MINLEN:50
		"""

rule get_18:
        output:
                'chr18.fa.gz'
        shell:
                """
                curl 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz' -o chr18.fa.gz
                """

rule df_18:
        input:
                'chr18.fa.gz'
        output:
                'chr18.fa'
        shell:
                '''
                gunzip {input}
                '''

rule get_18annot:
        output:
                'chr18.gtf'
        shell:
                """
                curl 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/GRCh37_mapping/gencode.v24lift37.basic.annotation.gtf.gz' | gunzip - -c > {output}
                """
                
rule index_18:
        input:
                'chr18.fa',
                'chr18.gtf'
        output:
                'chr18.idx/SA'
                
        shell:
                """
                STAR --runMode genomeGenerate \
                --genomeDir chr18.idx \
                --genomeFastaFiles {input[0]} \
                --sjdbGTFfile {input[1]}
                """
rule map:
	input:
		'trimmed/{seq}.R1P.fastq',
		'trimmed/{seq}.R2P.fastq',
                'chr18.idx/SA'
	output:
		'mapped/{seq}.Aligned.sortedByCoord.out.bam'
	shell:
		"""
		STAR --outFilterMultimapNmax 1\
		--genomeDir `dirname {input[2]}` \
		--outSAMattributes All --outSAMtype BAM SortedByCoordinate \
		--outFileNamePrefix mapped/{wildcards.seq}. \
		--readFilesIn {input[0]} {input[1]}
		"""

rule index_mappings:
    input:
        'mapped/{seq}.Aligned.sortedByCoord.out.bam'
    output:
        'mapped/{seq}.Aligned.sortedByCoord.out.bam.bai'
    shell:
        """
        samtools index {input}
        """
rule count_features:
    input:
        'chr18.gtf',
        bams
    output:
        'fcounts.txt'
    shell:
        """
        featureCounts -p -t exon -g gene_id \
                -o {output} -a {input}
        """

rule tabulate:
    input:
        'chr18.gtf',
        'fcounts.txt'
    output:
        'counts.txt'
    shell:
        """
        export LC_ALL=C
        join \
            <(perl -ne 'print "$1 $2\n" if /gene_id \"(.*?)\".*gene_name \"(.*?)\"/' \
            {input[0]} | sort | uniq | sort) \
            <(sort {input[1]}) \
            | awk '{{print $2 " " $7 " " $8 " " $9 " " $10 " " $11 " " $12}}' \
            > {output}
        """
