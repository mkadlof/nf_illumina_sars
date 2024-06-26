process freeBayes {
    tag "freeBayes:${sampleId}"
    publishDir "${params.results_dir}/${sampleId}/freebayes", mode: 'copy'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple val(sampleId2), path(ref_genome)

    output:
    tuple val(sampleId), path('freebayes.fa')

    script:
    """
    freebayes --limit-coverage ${params.max_depth} \
              --min-coverage ${params.min_cov} \
              --min-mapping-quality ${params.min_mapq} \
              --min-base-quality ${params.quality_snp} \
              --use-mapping-quality \
              --fasta-reference ${ref_genome} \
              --ploidy 1 \
              ${bam} > detected_variants_freebayes.vcf

    cat detected_variants_freebayes.vcf | \
        bcftools norm --check-ref w \
                      --rm-dup all \
                      --fasta-ref ${ref_genome} | \
                          bcftools norm --check-ref w \
                                        --multiallelics -indels \
                                        --fasta-ref ${ref_genome} > detected_variants_freebayes_fix.vcf
    
    qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`
    
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >= ${params.min_cov} & (SAF  + SAR)/(SRF + SRR + SAF + SAR) > ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > detected_variants_freebayes_fix_high.vcf
    
    bgzip --force detected_variants_freebayes_fix_high.vcf
    tabix detected_variants_freebayes_fix_high.vcf.gz
    
    bcftools filter --include "QUAL >= \${qual} & INFO/DP >=  ${params.min_cov}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) >= ${params.lower_ambig}  & (SAF  + SAR)/(SRF + SRR + SAF + SAR) <= ${params.upper_ambig} " \
             detected_variants_freebayes_fix.vcf > tmp_low.vcf
    introduce_amb_2_vcf.py tmp_low.vcf \
           detected_variants_freebayes_fix_ambig.vcf

    bgzip --force detected_variants_freebayes_fix_ambig.vcf
    tabix detected_variants_freebayes_fix_ambig.vcf.gz

    bcftools concat detected_variants_freebayes_fix_high.vcf.gz \
                    detected_variants_freebayes_fix_ambig.vcf.gz | \
                        bcftools sort --output-type z > detected_variants_freebayes_final.vcf.gz
    tabix detected_variants_freebayes_final.vcf.gz

    cat ${ref_genome} | \
            bcftools consensus --mark-del X --samples - \
                               detected_variants_freebayes_final.vcf.gz > freebayes.fa
    """
}
