process detect_subtype {
    tag "detect_subtype:${sampleId}"
    maxForks 5

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("subtype_counts_each_segment.txt"), path("subtype_scores_each_segment.txt"), env(REF_GENOME_ID)

    script:
    """
    run_bwa() {
        for GENOME in "\${@}"; do
            bwa mem -t ${params.threads} -T 30 /home/data/infl/genomes/\${GENOME}/\${GENOME}.fasta ${reads[0]} ${reads[1]} | \
                samtools view -@ ${params.threads} -Sb -f 3 -F 2048 - | \
                samtools sort -@ ${params.threads} -o \${GENOME}.bam -
            samtools index \${GENOME}.bam
        done
    }

    find_max() {
        local GENOMES=("\${!1}")
        local VALUES=("\${!2}")
        local max=0
        local max_index=0
        local i=0
        for value in \${VALUES[@]}
        do
            echo \${GENOMES[\${i}]} \${value} >> subtype_scores.txt
            if [ \${value} -gt \${max} ]; then
                max=\${value}
                max_index=\${i}
            fi
            ((i++))
        done
        echo \${GENOMES[\${max_index}]}
    }

    variant="${params.variant}"

    KNOWN_VARIANTS='H1N1 H3N2 H4N6 H5N2 H5N1 H5N6 H5N8 H6N1 H7N9 H9N2 Yamagata Victoria UNK'
    ALL_GENOMES=(`ls /home/data/infl/genomes`)
    ALL_SEGMENTS=(PB2 PB1 PA HA NP NA MP NS)

    echo -e "id \${ALL_SEGMENTS[@]}" | tr " " "\t" >> subtype_counts_each_segment.txt
    echo -e "id \${ALL_SEGMENTS[@]}" | tr " " "\t" >> subtype_scores_each_segment.txt

    if [ \${variant} == 'UNK' ] ;then
        # We map our sequences to each of the references and count the number of reads mapping to
        # the HA and NA segments. The reference with the most reads wins.

        run_bwa \${ALL_GENOMES[@]}

        for GENOME in \${ALL_GENOMES[@]}; do
            alignment_score+=(`get_alignment_score.py \${GENOME}.bam chr6_NA,chr4_HA`)
            echo -e "\${GENOME}\t`get_alignment_score_all_segments.py \${GENOME}.bam`" >> subtype_scores_each_segment.txt
            MEAN_COVERAGE=()
            for SEGMENT in \${ALL_SEGMENTS[@]}; do
                MEAN_COVERAGE+=(`bedtools genomecov -ibam \${GENOME}.bam  -bga | \
                                    grep -w chr[0-9]_\${SEGMENT} | \
                                    awk 'BEGIN {segmentu_length=0; read_count=0} {segment_length+=(\$3-\$2); read_count+=\$4 * (\$3-\$2)} END {print read_count/segment_length}'`)
            done
            echo -e "\${GENOME} \${MEAN_COVERAGE[@]}" | tr " " "\t" >> subtype_counts_each_segment.txt
        done

        result=`find_max ALL_GENOMES[@] alignment_score[@]`
        echo \${result}
        result_mini=`echo \${result} | cut -d "_" -f1`

        PRIMERS="data/infl/primers/\${result}/\${result}_primers.bed"
        REF_GENOME_FASTA="data/infl/genomes/\${result}/\${result}.fasta"
    else
        SAMPLE_GENOMES=()

        for genome in "\${ALL_GENOMES[@]}"; do
            if [[ "\${genome}" == *"\${variant}"* ]]; then
                    SAMPLE_GENOMES+=(\${genome})
            fi
        done

        run_bwa \${SAMPLE_GENOMES[@]}

        alignment_score=()
        for GENOME in \${SAMPLE_GENOMES[@]}; do
            alignment_score+=(`get_alignment_score.py \${GENOME}.bam chr6_NA,chr4_HA`)
            echo -e "\$GENOME\t`get_alignment_score_all_segments.py \${GENOME}.bam`" >> subtype_scores_each_segment.txt
            MEAN_COVERAGE=()
            for SEGMENT in \${WSZYSTKIE_SEGMENTY[@]}; do
                MEAN_COVERAGE+=(`bedtools genomecov -ibam  \${GENOME}.bam  -bga | grep -w chr[0-9]_\${SEGMENT} | awk 'BEGIN {segment_length=0; read_count=0} {segment_length+=(\$3-\$2); read_count+=\$4 * (\$3-\$2)} END {print read_count/segment_length}'`)
            done
            echo -e "\${GENOME} \${MEAN_COVERAGE[@]}" | tr " " "\t" >> subtype_counts_each_segment.txt
        done

        result=`find_max SAMPLE_GENOMES[@] alignment_score[@]`
        result_mini=`echo \${result} | cut -d "_" -f1`
        PRIMERS="data/infl/primers/\${result}/\${result}_primers.bed"
        REF_GENOME_FASTA="data/infl/genomes/\${result}/\${result}.fasta"
    fi
    REF_GENOME_ID=\${result}
    rm *bam*
    """
}