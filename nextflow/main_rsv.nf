#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// include modules
include {printHelp} from "${params.nfpath}/nextflow/modules/help.nf"
include {makeFastqSearchPath} from "${params.nfpath}/nextflow/modules/util.nf"

// import subworkflows
include {denovo} from "${params.nfpath}/nextflow/workflows/workflow_rsv.nf"

def raiseError ( value ) {
    sleep(2000)
    println(value)
    System.exit(1)
}

if (params.help) {
    printHelp()
    exit 0
}

// main workflow
workflow {
    if ( params.directory ) {
        fastqSearchPath = makeFastqSearchPath( params.illuminaPrefixes, params.illuminaSuffixes, params.fastq_exts )

        Channel.fromFilePairs( fastqSearchPath, flat: true)
            .filter{ !( it[0] =~ /Undetermined/ ) }
            .set{ ch_filePairs }
    }
    else if( params.fastq ) {
        Channel.fromPath(params.fastq).toSortedList().flatten().buffer(size: 2).set{ ch_filePairs }
        ch_filePairs.map{ it -> [ (it[1] =~ /(.+)\/([^\/]+)_R([^\._]+)(.+)/)[0][2], it[0], it[1] ] }.set{ ch_filePairs }
    }
    else {
        println("Please provide the path of the directory containing your fastq pairs or provide each file's path in the json configuration file.")
        System.exit(1)
    }

    ch_filePairs.map{ (it[1] =~ /(.+)\/([^\/]+)_R([^\._]+)(.+)/)[0][2] != (it[2] =~ /(.+)\/([^\/]+)_R([^\._]+)(.+)/)[0][2] ? raiseError("Fastq file pairs are malformed. Error raised for " + (it[1] =~ /(.+)\/([^\/]+)_R([^\._]+)(.+)/)[0][2] + "/" + (it[2] =~ /(.+)\/([^\/]+)_R([^\._]+)(.+)/)[0][2] + "."): '' }

    main:
    if ( params.workflow == 'denovo') {
        denovo(ch_filePairs)
    }
    
}
