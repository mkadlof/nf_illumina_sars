process varScan {
    publishDir "results/${sampleId}", mode: 'symlink'

    input:
    tuple val(sampleId), path(bam), path(bai)
    tuple path(reference_fasta), path(reference_fai)

    output:
    tuple val(sampleId), path('detected_variants_varscan_final.vcf.gz'), path('detected_variants_varscan_final.vcf.gz.tbi')
    tuple val(sampleId), path('varscan.fa')

    script:
    """
    samtools mpileup --max-depth ${params.max_depth} \
                 --fasta-ref ${reference_fasta} \
                 --min-BQ ${params.quality_snp} \
                 ${bam} >> ${bam}.mpileup

    varscan_qual=`echo "${params.quality_snp} - 1" | bc -l`
    java -jar /opt/varscan/VarScan.v2.4.6.jar pileup2cns ${bam}.mpileup \
                                                         --min-avg-qual \${varscan_qual} \
                                                         --p-value ${params.pval} \
                                                         --min-var-freq ${params.lower_ambig} \
                                                         --min-coverage ${params.min_cov} \
                                                         --variants \
                                                         --min-reads2 0 > detected_variants_varscan.txt

    parse_vcf_output_final.py detected_variants_varscan.txt ${params.upper_ambig} ${params.pval}

    bgzip --force detected_variants_varscan.vcf
    tabix detected_variants_varscan.vcf.gz

    qual=`echo ${params.pval} | awk '{print int(10*-log(\$1)/log(10))}'`

    bcftools norm --check-ref w \
                  --rm-dup all \
                  --fasta-ref ${reference_fasta}\
                   detected_variants_varscan.vcf.gz | \
                       bcftools norm --check-ref w \
                                     --multiallelics -indels \
                                     --fasta-ref ${reference_fasta} | \
                                           bcftools filter \
                                                    --include "QUAL >= \${qual} && AF >= ${params.lower_ambig} && DP >= ${params.min_cov}" > detected_variants_varscan_final.vcf

    bgzip --force detected_variants_varscan_final.vcf
    tabix detected_variants_varscan_final.vcf.gz

    cat ${reference_fasta} | bcftools consensus --samples - detected_variants_varscan_final.vcf.gz > varscan.fa
    """
}