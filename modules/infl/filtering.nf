process filtering {
    tag "filtering:${sampleId}"

    input:
    tuple val(sampleId), path(bam), path(bai), path(ref_genome_fasta), path(primers)

    output:
    tuple val(sampleId), path('to_clip_sorted.bam'), path('to_clip_sorted.bam.bai')

    script:
    """
    # Custom filtering for influenza
    # The script takes in sequence:
    # - BAM for filtering and downsampling
    # - Primer scheme
    # - Target coverage per segment (integer)
    # - Minimum read mapping quality (integer)
    # - Reference genome sequences in FASTA (all segments)

    simple_filter_illumina_INFL.py ${bam} ${primers} ${params.max_depth} ${params.min_mapq} ${params.length} ${ref_genome_fasta}
    """
}
