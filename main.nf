#!/usr/bin/env nextflow

// Edit nextflow.configuration!

data_location=config.data_location
large_core=config.large_core
small_core=config.small_core

// ** - Recurse through subdirectories
fq_set = Channel.fromPath(data_location + "**/*.fastq.gz")
                .map { n -> [ n.getName(), n ] }


project="PRJNA13758"
reference="WS255"
prefix="http://ftp.wormbase.org/pub/wormbase/releases/${reference}/species/c_elegans/${project}"


process fetch_reference {

    publishDir "output/", mode: 'copy', pattern: 'meta.pipeline.txt'
    
    output:
        file("geneset.gtf.gz") into geneset_gtf

        file("reference.fa.gz") into reference_hisat
        file("meta.pipeline.txt")

    """
        echo '${project}' > meta.pipeline.txt
        echo '${reference}' >> meta.pipeline.txt
        # Download gene set
        curl ${prefix}/c_elegans.${project}.${reference}.canonical_geneset.gtf.gz > geneset.gtf.gz

        # Download reference genome
        curl ${prefix}/c_elegans.${project}.${reference}.genomic.fa.gz > reference.fa.gz
    """
}
geneset_gtf.into { geneset_hisat; geneset_stringtie }


extract_exons_py = file("scripts/hisat2_extract_exons.py")
extract_splice_py = file("scripts/hisat2_extract_splice_sites.py")

process hisat2_indexing {

    input:
        file("geneset.gtf.gz") from geneset_hisat
        file("reference.fa.gz") from reference_hisat

    output:
        file("splice.ss") into splice_hisat
        file("exon.exon") into exon_hisat
        file("reference.fa.gz") into reference_build_hisat

    """
        gzcat geneset.gtf.gz | python ${extract_splice_py} - > splice.ss
        gzcat geneset.gtf.gz | python ${extract_exons_py} - > exon.exon
    """

}

process build_hisat_index {

    cpus small_core

    input:
        file("splice.ss") from splice_hisat
        file("exon.exon") from exon_hisat
        file("reference.fa.gz") from reference_build_hisat

    output:
        file "*.ht2" into hs2_indices

    """
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${small_core} --ss splice.ss --exon exon.exon reference.fa reference.hisat2_index
    """

}

process trimmomatic {

    cpus small_core

    tag { name }

    input:
        set val(name), file(reads) from fq_set

    output:
        file(name_out) into trimmed_reads

    script:
    name_out = name.replace('.fastq.gz', '_trim.fq.gz')

    """
        trimmomatic SE -phred33 -threads ${small_core} ${reads} ${name_out} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15 
    """
}


process align {

    cpus small_core

    tag { prefix }

    input:
        file reads from trimmed_reads
        file hs2_indices from hs2_indices.first()

    output:
        set val(sample_id), file("${prefix}.bam"), file("${prefix}.bam.bai") into hisat2_bams
        file "${prefix}.hisat2_log.txt" into alignment_logs

    script:
        index_base = hs2_indices[0].toString() - ~/.\d.ht2/
        prefix = reads[0].toString() - ~/(_trim)(\.fq\.gz)$/
        m = prefix =~ /\w+-([^_]+)_.*/
        sample_id = m[0][1]

    """
        hisat2 -p ${small_core} -x $index_base -U ${reads} -S ${prefix}.sam --rg-id "${prefix}" --rg "SM:${sample_id}" --rg "LB:${sample_id}" --rg "PL:ILLUMINA" 2> ${prefix}.hisat2_log.txt
        samtools view -bS ${prefix}.sam > ${prefix}.unsorted.bam
        samtools flagstat ${prefix}.unsorted.bam
        samtools sort -@ ${small_core} -o ${prefix}.bam ${prefix}.unsorted.bam
        samtools index -b ${prefix}.bam
    """
}

joint_bams = hisat2_bams.groupTuple()

process join_bams {

    publishDir "output/bam", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from joint_bams


    output:
        set val(sample_id), file("${sample_id}.bam"), file("${sample_id}.bam.bai") into hisat2_merged_bams

    """
        samtools merge -f ${sample_id}.unsorted.bam ${bam}
        samtools sort -@ ${small_core} -o ${sample_id}.bam ${sample_id}.unsorted.bam
        samtools flagstat ${sample_id}.bam
        samtools index -b ${sample_id}.bam

    """
}


process stringtie_counts {

    //storeDir '/storedir/expression'
    publishDir "output/expression", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        set val(sample_id), file(bam), file(bai) from hisat2_merged_bams
        file("geneset.gtf.gz") from geneset_stringtie.first()

    output:
        file("${sample_id}/*") into stringtie_exp

    """ 
        gzcat geneset.gtf.gz > geneset.gtf
        stringtie -p ${small_core} -G geneset.gtf -e -B -o ${sample_id}/${sample_id}_expressed.gtf ${bam}
    """
}


prepDE = file("scripts/prepDE.py")

process stringtie_table_counts {

    echo true

    publishDir "output/diffexp", mode: 'copy'

    cpus small_core

    tag { sample_id }

    input:
        val(sample_file) from stringtie_exp.toSortedList()

    output:
        file ("gene_count_matrix.csv") into gene_count_matrix
        file ("transcript_count_matrix.csv") into transcript_count_matrix

    """
        for i in ${sample_file.flatten().join(" ")}; do
            bn=`basename \${i}`
            full_path=`dirname \${i}`
            sample_name=\${full_path##*/}
            echo "\${sample_name} \${i}"
            mkdir -p expression/\${sample_name}
            ln -s \${i} expression/\${sample_name}/\${bn}
        done;
        python ${prepDE} -i expression -l 50 -g gene_count_matrix.csv -t transcript_count_matrix.csv

    """
}



