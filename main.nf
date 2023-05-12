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

// run trim galore

process trimgalore {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus params.trimgaloreJobCpus

    tag "Running trim galore on the samples"
    publishDir path:{params.outputDir},mode: 'symlink'

    input:
    file sampleID from concatfastq1_ch

    output:
    file("${sampleID}") into trimmedfastq1_ch
    file("${sampleID}") into trimmedfastq2_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    trim_galore --cores ${params.trimgaloreJobCpus} -o ${sampleID} --paired ${sampleID}/${sampleID}_1.fastq.gz ${sampleID}/${sampleID}_2.fastq.gz

    mv ${sampleID}/${sampleID}_1_val_1.fq.gz ${sampleID}/${sampleID}.trimmed_1.fastq.gz

    mv ${sampleID}/${sampleID}_2_val_2.fq.gz ${sampleID}/${sampleID}.trimmed_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """
}

// run fastqc on raw fastq

process raw_fastqc {
    container 'munozcriollojj/nf-pipeline-test:latest'

    cpus 1

    tag "Running fastQC on raw fastq"
    publishDir path:{params.outputDir},mode: 'symlink'

    input:
    file sampleID from concatfastq2_ch

    output:
    file("${sampleID}") into fastqonraw_ch

    script:
    """ 
    sleep ${params.sleepTimeStart}

    fastqc -o ${sampleID} ${sampleID}/${sampleID}_1.fastq.gz

    fastqc -o ${sampleID} ${sampleID}/${sampleID}_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """
}

// run fastqc on trimmed fastq

process trimmed_fastqc {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus 1

    tag "Running fastQC on trimmed fastq"
    publishDir path:{params.outputDir},mode: 'symlink'

    input:
    file sampleID from trimmedfastq1_ch

    output:
    file("${sampleID}") into fastqontrimmed_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_1.fastq.gz

    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_2.fastq.gz

    sleep ${params.sleepTimeEnd}
    """
}

// run multiqc on all fastqc results

process multiqc {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus 1

    tag "Running multiqc on fastqc reports"
    publishDir path:{params.outputDir},mode: 'symlink'

    input:
    file dummy from fastqontrimmed_ch.collect()

    output:
    file("multiQC.html") into multiqc1_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${params.projectDir}/trace.txt ${params.projectDir}/work/ ${params.projectDir}/${params.outputDir} \$linkDir trimmed_fastqc

    multiqc \$linkDir -n multiQC

    rm -r \$linkDir

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
