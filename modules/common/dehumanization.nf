process dehumanization_illumina  {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path('forward_paired_nohuman.fastq.gz'), path('reverse_paired_nohuman.fastq.gz'), emit: to_pubdir
    tuple val(sampleId), path('dehumanization.json'), emit: json

    script:
    """
    if [ ${QC_status} == "nie" ]; then
      touch forward_paired_nohuman.fastq.gz 
      touch reverse_paired_nohuman.fastq.gz
      touch dehumanization.json
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads[0]} lista_id_nohuman.txt >> forward_paired_nohuman.fastq
      seqtk subseq ${reads[1]} lista_id_nohuman.txt >> reverse_paired_nohuman.fastq

      gzip forward_paired_nohuman.fastq
      gzip reverse_paired_nohuman.fastq
      # For now dummy json
      touch dehumanization.json
    fi
    """
}

process dehumanization_nanopore {
    tag "dehumanization:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy'

    input:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai'), val(QC_status), path(reads)

    output:
    tuple val(sampleId), path('forward_nohuman.fastq.gz'), emit: to_pubdir
    tuple val(sampleId), path('dehumanization.json'), emit: json
    script:
    """

    if [ ${QC_status} == "nie" ]; then
      touch forward_nohuman.fastq.gz
      touch dehumanization.json
    else
      samtools view mapped_reads.bam | cut -f1 | sort | uniq >> lista_id_nohuman.txt
      seqtk subseq ${reads} lista_id_nohuman.txt >> forward_nohuman.fastq

      gzip forward_nohuman.fastq
      # For now dummy json
      touch dehumanization.json
    fi
    """
}
