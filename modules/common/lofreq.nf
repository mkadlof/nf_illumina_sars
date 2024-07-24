process lofreq {
    tag "lofreq:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/lofreq", mode: 'copy'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId2), path(ref_genome_with_index)

    output:
    tuple val(sampleId), path('lofreq.fa')

    script:
    """
    # If index file is missing in the input, then create it right now.
    index_file="${ref_genome_with_index[0]}.fai"
    if [ ! -f "\$index_file" ]; then
        samtools faidx "${ref_genome_with_index[0]}"
    fi

    lofreq call-parallel --pp-threads ${params.threads} \
                         --ref ${ref_genome_with_index[0]} \
                         --max-depth ${params.max_depth} \
                         --min-cov ${params.min_cov} \
                         --call-indels \
                         --out detected_variants_lofreq.vcf \
                         ${bam}
    
    cat detected_variants_lofreq.vcf | \
            bcftools norm --check-ref w \
                          --rm-dup all \
                          --fasta-ref ${ref_genome_with_index[0]} | \
                              bcftools norm --check-ref w \
                                            --multiallelics -indels \
                                            --fasta-ref ${ref_genome_with_index[0]} > detected_variants_lofreq_fix.vcf

    qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`
    cat detected_variants_lofreq_fix.vcf | \
        awk -v qual=\${qual} \
            -v min_cov=${params.min_cov} \
            -v upper_ambig=${params.upper_ambig} \
            -v lower_ambig=${params.lower_ambig} \
            '{if(substr(\$1, 1, 1) == "#") print \$0; else split(\$8,a, ";"); split(a[4], b, "="); split(b[2], c, ","); split(a[1], DP, "="); split(a[2], AF, "="); if( ((c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  > upper_ambig &&  \$6 >= qual  && DP[2] >= min_cov) || (a[5] == "INDEL" && AF[2] >= lower_ambig && \$6 >= qual  && DP[2] >= min_cov)) print \$0}' > detected_variants_lofreq_fix_high.vcf
    
    bgzip --force detected_variants_lofreq_fix_high.vcf
    tabix detected_variants_lofreq_fix_high.vcf.gz
    
    cat detected_variants_lofreq_fix.vcf | \
        awk -v qual=\${qual} \
            -v min_cov=${params.min_cov} \
            -v upper_ambig=${params.upper_ambig} \
            -v lower_ambig=${params.lower_ambig} \
            '{if(substr(\$1, 1, 1) == "#") print \$0; else split(\$8,a, ";"); split(a[4], b, "="); split(b[2], c, ","); split(a[1], DP, "="); if( (c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  <= upper_ambig &&  \$6 >= qual  && DP[2] >= min_cov && (c[3] + c[4]) / (c[3] + c[4] + c[1] + c[2])  >= lower_ambig && (a[5] != "INDEL")) print \$0}' > tmp_lofreq.vcf
    
    introduce_amb_2_vcf.py tmp_lofreq.vcf\
                                                 detected_variants_lofreq_fix_ambig.vcf
    bgzip --force detected_variants_lofreq_fix_ambig.vcf
    tabix detected_variants_lofreq_fix_ambig.vcf.gz
    
    bcftools concat detected_variants_lofreq_fix_high.vcf.gz \
                    detected_variants_lofreq_fix_ambig.vcf.gz | \
                        bcftools sort --output-type z > detected_variants_lofreq_final.vcf.gz
    tabix detected_variants_lofreq_final.vcf.gz
    
    cat ${ref_genome_with_index[0]} | \
        bcftools consensus --mark-del X --samples - \
                           detected_variants_lofreq_final.vcf.gz > lofreq.fa
    """
}




