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

    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir trimmed_fastqc

    multiqc \$linkDir -n multiQC

    rm -r \$linkDir

    sleep ${params.sleepTimeEnd}
    """
}

// get sequence length for STAR parameter

process sequencelength {
    container 'munozcriollojj/nf-pipeline-test:latest'

    tag "Calculating sequence length for STAR genome index"
    publishDir path:{params.resourcesDir},mode: 'symlink'

    output:
    file("seqlen.txt") into sequencelength_ch

    script:
    """
    bash ${params.projectDir}/${params.srcDir}/calculate_fastq_sequence_length.sh ${params.dataDir} > seqlen.txt
    """
}

// copy genome assembly file into resources

process copy_genome {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus 1

    tag "Copying genome file into resources/"
    publishDir path:{params.resourcesDir},mode: 'symlink'

    output:
    file("${params.genomeName}") into genome_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    cp ${launchDir}/resources/assembly/${params.genomeName} .

    sleep ${params.sleepTimeEnd}
    """
}

// copy GTF file into resources

process copy_gtf {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus 1

    tag "Copying gtf file into resources/"
    publishDir path:{params.resourcesDir},mode: 'symlink'

    output:
    file("*") into gtf2_ch
    file("*") into gtf3_ch
    file("*") into gtf4_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    cp ${launchDir}/resources/gtf/${params.gtfName} .

    sleep ${params.sleepTimeEnd}
    """
}

// create STAR genome index

process genome_index {
    container 'munozcriollojj/nf-pipeline-test:latest'
    cpus params.genomeIndexJobCpus

    tag "Indexing genome using STAR"
    publishDir path:{params.resourcesDir},mode: 'symlink'

    input:
    file gtf from gtf2_ch
    file genome from genome_ch
    file seqlen from sequencelength_ch

    output:
    file "genome_index" into genomeindex_ch

    script:
    """
    sleep ${params.sleepTimeStart}

    mkdir genome_index

    STAR --runThreadN ${params.genomeIndexJobCpus} --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf} --sjdbOverhang `cat ${seqlen}`

    sleep ${params.sleepTimeEnd}
    """
}
