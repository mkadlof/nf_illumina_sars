process freeBayes {

    input:
    tuple val(sampleId), path(bam), path(bai), path(lowcoverage_masked_fa)
    path ref_genome

    output:
    tuple val(sampleId), path('detected_variants_freebayes_final.vcf.gz')
    tuple val(sampleId), path('freebayes_masked.fa')

    script:
    """
    freebayes --limit-coverage ${params.max_depth} \
              --min-coverage ${params.min_cov} \
              --min-mapping-quality 20 \
              --min-base-quality ${params.quality_SNP} \
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
            bcftools consensus --samples - \
                               detected_variants_freebayes_final.vcf.gz > freebayes.fa
    cat ${lowcoverage_masked_fa} \
        freebayes.fa > tmp_freebayes.fa

    mafft --auto --inputorder --quiet tmp_freebayes.fa > tmp_freebayes_aln.fa

    get_N.py tmp_freebayes_aln.fa
    
    mv output_freebayes_masked.fa freebayes_masked.fa
    """
}