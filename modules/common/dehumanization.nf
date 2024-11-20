process dehumanization_illumina  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'
    container  = params.main_image

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('list_of_dehumanzed_fastas.txt'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_forward_paired_nohuman.fastq.gz ${sampleId}_reverse_paired_nohuman.fastq.gz
    else
        samtools view mapped_reads.bam | cut -f1 | sort | uniq > lista_id_nohuman.txt
        seqtk subseq ${reads[0]} lista_id_nohuman.txt | gzip > ${sampleId}_forward_paired_nohuman.fastq.gz
        seqtk subseq ${reads[1]} lista_id_nohuman.txt | gzip > ${sampleId}_reverse_paired_nohuman.fastq.gz
    fi

    ls *paired_nohuman* | tr " " "\n" >> list_of_dehumanzed_fastas.txt
    """
}

process dehumanization_nanopore {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'
    container  = params.main_image
    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path("*nohuman.fastq.gz"), emit: to_pubdir
    tuple val(sampleId), path('list_of_dehumanzed_fastas.txt'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch ${sampleId}_nohuman.fastq.gz
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads} lista_id_nohuman.txt | gzip >> ${sampleId}_nohuman.fastq.gz
    fi

    ls *nohuman* >> list_of_dehumanzed_fastas.txt
    """
}
