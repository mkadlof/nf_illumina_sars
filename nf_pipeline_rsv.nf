// Programs parameters can be modified by a shell wrapper
params.memory = 2024
params.min_number_of_reads = 0
params.min_median_quality = 0
params.quality_initial = 5
params.length = 90
params.max_number_for_SV = 200000
params.max_depth = 600
params.min_cov = 20
params.mask = 20
params.quality_snp = 15
params.pval = 0.05
params.lower_ambig = 0.45
params.upper_ambig = 0.55
params.window_size = 50 // Window size in which we equalize the coverage
params.min_mapq = 30


// Specifying location of reads and primers and necessary databases, MUST be selected by a user
params.threads = 5
params.reads = "" // Must be provided by user
params.primers_id = "" // We replace the name of the available directory in the wrapper with V4, V4.1, etc. The user provides this.
params.pangolin_db_absolute_path_on_host = ""
params.nextclade_db_absolute_path_on_host = ""
params.kraken2_db_absolute_path_on_host = ""
params.freyja_db_absolute_path_on_host = ""

// adapters, can be set  by a user but there is a default
params.adapters_id="TruSeq3-PE-2" // To podaje user ale jest default

// output dir, by default "results" directory
params.results_dir = "./results"

// Directory with modules, maybe move the .nf file to the main script to remove one settable parameter ? MUST be indicated by user
params.modules = "/absolute/path/to/directory/with/modules"

adapters="/home/data/common/adapters/${params.adapters_id}.fa"

params.ref_genome_id="MN908947.3"

include { kraken2 } from "${params.modules}/common/kraken2.nf"
include { bwa } from "${params.modules}/common/bwa.nf"
include { dehumanization } from "${params.modules}/common/dehumanization.nf"
include { fastqc as fastqc_1 } from "${params.modules}/common/fastqc.nf"
include { fastqc as fastqc_2 } from "${params.modules}/common/fastqc.nf"
include { trimmomatic } from "${params.modules}/common/trimmomatic.nf"
include { filtering } from "${params.modules}/sarscov2/filtering.nf"
include { masking } from "${params.modules}/common/masking.nf"
include { merging } from "${params.modules}/sarscov2/merging.nf"
include { picard } from "${params.modules}/common/picard.nf"
include { manta } from "${params.modules}/common/manta.nf"
include { indelQual } from "${params.modules}/common/indelQual.nf"
include { wgsMetrics } from "${params.modules}/common/wgsMetrics.nf"
include { lowCov } from "${params.modules}/common/lowCov.nf"
include { varScan } from "${params.modules}/common/varscan.nf"
include { freeBayes } from "${params.modules}/common/freeBayes.nf"
include { lofreq } from "${params.modules}/common/lofreq.nf"
include { consensus } from "${params.modules}/common/consensus.nf"
include { vcf_for_fasta } from "${params.modules}/sarscov2/vcf_for_fasta.nf"
include { snpEff } from "${params.modules}/sarscov2/snpEff.nf"
include { simpleStats } from "${params.modules}/sarscov2/simpleStats.nf"
include { nextclade } from "${params.modules}/sarscov2/nextclade.nf"
include { pangolin } from "${params.modules}/sarscov2/pangolin.nf"
include { modeller } from "${params.modules}/sarscov2/modeller.nf"
include { json_aggregator } from "${params.modules}/common/json_aggregator.nf"

// Coinfection line
include { freyja } from "${params.modules}/sarscov2/freyja.nf"
include { coinfection_ivar } from "${params.modules}/sarscov2/coinfection_ivar.nf"
include { coinfection_varscan } from "${params.modules}/sarscov2/coinfection_varscan.nf"
include { coinfection_analysis } from "${params.modules}/sarscov2/coinfection_analysis.nf"

workflow{
    // Channels
    reads = Channel.fromFilePairs(params.reads)
    ref_genome = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta", size: -1).collect().sort()
    ref_genome_with_index = Channel.fromFilePairs("${projectDir}/data/sarscov2/genome/*sarscov2.fasta{,.amb,.ann,.bwt,.fai,.pac,.sa}", size: -1).collect().sort()
    primers_and_pairs = Channel.fromFilePairs("${projectDir}/data/sarscov2/primers/${params.primers_id}/*{nCoV-2019.scheme.bed,pairs.tsv}").collect()
    coinfections = Channel.fromPath("${projectDir}/data/sarscov2/coinfections/").first()
    vcf_template = Channel.fromPath("${projectDir}/data/sarscov2/vcf_template/").first()
    modeller_data = Channel.fromPath("${projectDir}/data/sarscov2/modeller/").first()

    // The following construct is needed to add the sampleId to the channels
    // with constant data (like reference genome or primers). This is done to
    // generalize cases where the reference genome may not be identical for
    // each sample (as is the case in the influenza pipeline, for example).
    ref_genome = ref_genome
    .combine(reads).map { _, genome, sampleId, reads ->
        return [sampleId, genome]
    }
    ref_genome_with_index = ref_genome_with_index
    .combine(reads).map { _, genome_and_index_files, sampleId, reads ->
        return [sampleId, genome_and_index_files]
    }
    primers_and_pairs = primers_and_pairs
    .combine(reads).map { _, primers_and_pairs, sampleId, reads ->
        return [sampleId, primers_and_pairs]
    }

    primers = primers_and_pairs.map{ sampleId, files -> [sampleId, files[0]] }

    // The following two variables are used exclusively to include pipeline version information in the resulting output.json file.
    def repo_path = workflow.projectDir
    version = Channel.value("git -C ${repo_path} rev-parse HEAD".execute().text.trim().substring(0, 7))
    pathogen = Channel.value('sars2')

    // Processes
    fastqc_1(reads, "initialfastq")
    c1 = reads.join(fastqc_1.out[2])
    kraken2(c1)
    trimmomatic(reads, adapters)
    fastqc_2(trimmomatic.out[0], "aftertrimmomatic")
    bwa(trimmomatic.out[0].join(ref_genome_with_index))
    c2 = bwa.out.join(trimmomatic.out[1])
    dehumanization(c2)
    c3 = bwa.out.join(primers)
    filtering(c3)
    c4 = filtering.out[0].join(primers_and_pairs)
    masking(c4)
    c5 = filtering.out[1].join(masking.out)
    merging(c5)
    picard(bwa.out)
    indelQual(merging.out.join(ref_genome))
    lowCov(indelQual.out.join(ref_genome))
    varScan(indelQual.out.join(ref_genome))
    freeBayes(indelQual.out.join(ref_genome))
    lofreq(indelQual.out.join(ref_genome_with_index))
    wgsMetrics(indelQual.out.join(ref_genome))
    c6 = lowCov.out[1].join(varScan.out).join(freeBayes.out).join(lofreq.out)
    consensus(c6)
    c7 = consensus.out[0].join(ref_genome)
    vcf_for_fasta(c7, vcf_template)
    manta(picard.out.join(consensus.out[0]))
    nextclade(manta.out)
    modeller(nextclade.out[1], modeller_data)
    pangolin(manta.out)
    c8 = vcf_for_fasta.out.join(indelQual.out).join(ref_genome)
    snpEff(c8)
    c9 = manta.out.join(wgsMetrics.out[0]).join(primers)
    simpleStats(c9)

    // Coinfection line
    c10 = bwa.out.join(ref_genome).join(primers)
    coinfection_ivar(c10)
    c11 = coinfection_ivar.out[0].join(ref_genome)
    freyja(c11)
    coinfection_varscan(coinfection_ivar.out[1])
    coinfection_analysis(coinfection_varscan.out, coinfections)


    json_aggregator(pathogen, version, wgsMetrics.out[1])
}
