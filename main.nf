#!/usr/bin/env nextflow

// Edit nextflow.configuration!

data_location=config.data_location
large_core=config.large_core
small_core=config.small_core

// ** - Recurse through subdirectories
fq_set = Channel.fromPath(data_location + "**/*.fastq.gz")


project="PRJNA13758"
reference="WS255"
prefix="ftp://ftp.wormbase.org/pub/wormbase/releases/${reference}/species/c_elegans/${project}/"


process fetch_reference {

    publishDir "output/", mode: 'copy', pattern: 'meta.pipeline.txt'
    
    output:
        file("geneset.gtf.gz") into geneset_hisat
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
        file("ref.*") into hisat_index

    """
        zcat reference.fa.gz > reference.fa
        hisat2-build -p ${small_core} --ss splice.ss --exon exon.exon reference.fa ref.hisat2_index
    """

}

hisat_index.println()
