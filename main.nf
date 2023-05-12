#!/usr/bin/env nextflow

 println """\
===================================
 N F 0 0 3
nextflow run main.nf -with-trace
===================================
General Information:
--------------------
  Profile(s):                 ${workflow.profile}
Input Parameters:
-----------------
  Base Directory:             ${params.projectDir}
  FastQC path:                ${params.dataDir}

"""

Channel.fromPath("${params.dataDir}/targets.csv")
       .splitCsv(header:true)
       .map{ row-> tuple(row.analysisID) }
       .into { sampleNames }

// concatenate and rename the input fastq files

process concatentate_and_rename_fastq {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus 1

    tag "Concatenating and renaming input fastq files"
    publishDir path:{params.outputDir},mode: 'symlink'

    input:
    set sampleID, file(read) from sampleNames

    output:
    file("${sampleID}") into concatfastq1_ch
    file("${sampleID}") into concatfastq2_ch

    script:
    """
    ${params.projectDir}/${params.srcDir}/concatenate_and_rename_fastq_input.pl ${params.dataDir} ${params.dataDir}/targets.csv ${sampleID}

    sleep ${params.sleepTimeEnd}

    """
}



// nextflow.enable.dsl=2
// params.str = 'Hello world!'
// 
// process splitLetters {
//   container 'nextflow/nextflow:21.10.6'
//   output:
//     path 'chunk_*'
// 
//   """
//   printf '${params.str}' | split -b 6 - chunk_
//   """
// }
// 
// process convertToUpper {
//   container 'nextflow/nextflow:21.10.6'
//   input:
//     path x
//   output:
//     stdout
// 
//   """
//   cat $x | tr '[a-z]' '[A-Z]'
//   """
// }
// 
// workflow {
//   splitLetters | flatten | convertToUpper | view { it.trim() }
// }
