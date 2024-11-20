process modeller {
    tag "modeller:${sampleId}"
    container  = params.main_image 
    publishDir "${params.results_dir}/${sampleId}", mode: 'copy', pattern: "${sampleId}_spike.pdb"

    input:
    tuple val(sampleId), path(target_fasta), val(QC_status)

    output:
    tuple val(sampleId), path('alignment.pir'), path("${sampleId}_spike.pdb"), emit: to_pubdir
    tuple val(sampleId), path('list_of_pdbs.txt'), emit: json
    script:
    """
    modeller_data="/home/data/sarscov2/modeller/"
    if [[ ${QC_status} == "nie" || ${params.species} != "SARS-CoV-2" ]]; then
        touch alignment.pir
        touch ${sampleId}_spike.pdb
    else
        cp \${modeller_data}/7dwz.pdb .
        modpy.sh modeller_create_alignment.py ${target_fasta}
        modpy.sh modeller_build_model.py alignment.pir
        cp target.B99990001.pdb ${sampleId}_spike.pdb
    fi

    ls -l *pdb >> list_of_pdbs.txt
    """
}
