/**
Nextflow FRASER workflow

**/
nextflow.enable.dsl=2


process FraserCount {
    time '50h'
    memory '80 GB'
    cpus 6

    input:
        path frasercountR
        path samplesheet
        path bamFiles
        path baiFiles
    output:
        path "fraser_output"

    script: 
        """
        mkdir "fraser_output"
        Rscript "${frasercountR}" "${samplesheet}" "fraser_output"
        """
}


process MergeCounts {
    time '10h'
    memory '80 GB'
    cpus 1

    input:
        path ext_counts
        path mergescriptR
        path fraser_output
        val ext_amount_fraser
        

    output:

        path fraser_output

    script: 
        """
        Rscript ${mergescriptR} "${fraser_output}" "${ext_counts}" "${ext_amount_fraser}"
        """

}

process Fraser {
    time '16h'
    memory '80 GB'
    cpus 6

    publishDir "$params.output/fraser", mode: 'copy'

    input:
        path samplesheet
        path output
        path fraserR

    output:

        path "*.tsv"
        path output

    script: 
        """
        Rscript "${fraserR}" "${samplesheet}" "${output}" "${params.genomeReferenceBuild}"
        """
}