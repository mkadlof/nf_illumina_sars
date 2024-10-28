process snpEff {
    tag "snpEff:${sampleId}"
    container  = params.main_image
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "detected_variants_consensus_annotated.txt"

    input:
    tuple val(sampleId), path(consensus_vcf_gz), path(consensus_vcf_gz_tbi), path('forvariants.bam'), path('forvariants.bam.bai'), val(QC_status),  path(ref_genome)
    output:
    tuple val(sampleId), path('detected_variants_consensus_annotated.vcf.gz'), path('detected_variants_consensus_annotated.txt')

    script:
    """
    ### We use Freebayes to calculate coverage and variant usage from the consensus VCF, but first,
    ### we generate the required BED file.
    zcat consensus.vcf.gz | grep MN | grep -v contig | awk '{print \$1, \$2-1, \$2}' | tr " " "\t" >> tmp.bed

    ### freebayes for predefined regions
    freebayes -t tmp.bed --fasta-reference ${ref_genome} --min-coverage 1 --min-mapping-quality 1 --min-base-quality 1 --ploidy 1 forvariants.bam >> freebayes_consesus_only.vcf

    ### running SNPeff
    java -jar /opt/snpEff/snpEff.jar ann -noStats ${params.ref_genome_id} \
         ${consensus_vcf_gz} > detected_variants_consensus_annotated.vcf

    bgzip --force detected_variants_consensus_annotated.vcf
    tabix detected_variants_consensus_annotated.vcf.gz
    
    # extracting infromations from snpEff output
    bcftools query --format '%POS | %REF%POS%ALT| %ANN \n' \
               detected_variants_consensus_annotated.vcf.gz | \
                  cut -d "|" -f1,2,4,6,13 | \
                  tr "|" "\t" | \
                  awk  'BEGIN {OFS = "\t"} {if ( \$3 == "upstream_gene_variant" || \$3 == "downstream_gene_variant") {gene="."; aa="."} else {gene=\$4; aa=\$5}; print \$1, gene, \$2, aa, \$3}' > part1.txt

   ### extracting DP and allele usage from freebayes output
   bcftools query --format "%POS %INFO/DP %INFO/SRF %INFO/SRR %INFO/SAF %INFO/SAR %POS \n" freebayes_consesus_only.vcf | awk '{print \$7, \$2, (\$6+\$5)/(\$3+\$4+\$5+\$6)}' >> part2.txt

   ### merging two files if freebayes dosn't have info regarding a variant from consensus.vcf we
   ### put "-" in DP and AF columns
   join -a 1 -1 1 -2 1 -o1.1,1.2,1.3,1.4,1.5,2.2,2.3  -e '-' part1.txt part2.txt | sort -nk1 | cut -d " " -f2- | tr " " "\t" >> detected_variants_consensus_annotated.txt
   """
}
