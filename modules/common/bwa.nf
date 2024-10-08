process bwa {
    tag "bwa:${sampleId}"
    maxForks 5

    input:
    tuple val(sampleId), path(reads), path(ref_genome_with_index)

    output:
    tuple val(sampleId), path('mapped_reads.bam'), path('mapped_reads.bam.bai')

    script:
    """
    bwa mem -t ${params.threads} -T 30 ${ref_genome_with_index[0]} ${reads[0]} ${reads[1]} | \
        samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
        samtools sort -@ ${params.threads} -o mapped_reads.bam -
    samtools index mapped_reads.bam
    """
}
