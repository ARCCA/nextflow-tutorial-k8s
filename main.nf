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
    ${params.projectDir}/${params.srcDir}/concatenate_and_rename_fastq_input.pl \
        ${params.dataDir} \
        ${params.dataDir}/targets.csv \
        ${sampleID}
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
    trim_galore \
        --cores ${params.trimgaloreJobCpus} \
        -o ${sampleID} \
        --paired ${sampleID}/${sampleID}_1.fastq.gz \
        ${sampleID}/${sampleID}_2.fastq.gz

    mv ${sampleID}/${sampleID}_1_val_1.fq.gz \
       ${sampleID}/${sampleID}.trimmed_1.fastq.gz

    mv ${sampleID}/${sampleID}_2_val_2.fq.gz \
       ${sampleID}/${sampleID}.trimmed_2.fastq.gz
    """
}

//// run fastqc on raw fastq
//
//process raw_fastqc {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//
//    cpus 1
//
//    tag "Running fastQC on raw fastq"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from concatfastq2_ch
//
//    output:
//    file("${sampleID}") into fastqonraw_ch
//
//    script:
//    """ 
//    sleep ${params.sleepTimeStart}
//
//    fastqc -o ${sampleID} ${sampleID}/${sampleID}_1.fastq.gz
//
//    fastqc -o ${sampleID} ${sampleID}/${sampleID}_2.fastq.gz
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// run fastqc on trimmed fastq
//
//process trimmed_fastqc {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Running fastQC on trimmed fastq"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from trimmedfastq1_ch
//
//    output:
//    file("${sampleID}") into fastqontrimmed_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_1.fastq.gz
//
//    fastqc -o ${sampleID} ${sampleID}/${sampleID}.trimmed_2.fastq.gz
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// run multiqc on all fastqc results
//
//process multiqc {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Running multiqc on fastqc reports"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file dummy from fastqontrimmed_ch.collect()
//
//    output:
//    file("multiQC.html") into multiqc1_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")
//
//    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir trimmed_fastqc
//
//    multiqc \$linkDir -n multiQC
//
//    rm -r \$linkDir
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// get sequence length for STAR parameter
//
//process sequencelength {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//
//    tag "Calculating sequence length for STAR genome index"
//    publishDir path:{params.resourcesDir},mode: 'symlink'
//
//    output:
//    file("seqlen.txt") into sequencelength_ch
//
//    script:
//    """
//    bash ${params.projectDir}/${params.srcDir}/calculate_fastq_sequence_length.sh ${params.dataDir} > seqlen.txt
//    """
//}
//
//// copy genome assembly file into resources
//
//process copy_genome {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Copying genome file into resources/"
//    publishDir path:{params.resourcesDir},mode: 'symlink'
//
//    output:
//    file("${params.genomeName}") into genome_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    cp ${launchDir}/resources/assembly/${params.genomeName} .
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// copy GTF file into resources
//
//process copy_gtf {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Copying gtf file into resources/"
//    publishDir path:{params.resourcesDir},mode: 'symlink'
//
//    output:
//    file("*") into gtf2_ch
//    file("*") into gtf3_ch
//    file("*") into gtf4_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    cp ${launchDir}/resources/gtf/${params.gtfName} .
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// create STAR genome index
//
//process genome_index {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus params.genomeIndexJobCpus
//
//    tag "Indexing genome using STAR"
//    publishDir path:{params.resourcesDir},mode: 'symlink'
//
//    input:
//    file gtf from gtf2_ch
//    file genome from genome_ch
//    file seqlen from sequencelength_ch
//
//    output:
//    file "genome_index" into genomeindex_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    mkdir genome_index
//
//    STAR --runThreadN ${params.genomeIndexJobCpus} --runMode genomeGenerate --genomeDir genome_index --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf} --sjdbOverhang `cat ${seqlen}`
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// run STAR mapping
//
//process star_mapping {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus params.starMappingJobCpus
//
//    tag "Mapping trimmed fastq with STAR"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from trimmedfastq2_ch
//    file genomeIndex from genomeindex_ch
//
//    output:
//    file("${sampleID}") into bam_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    STAR --readFilesCommand zcat --runThreadN ${params.starMappingJobCpus} --runMode alignReads --outSAMunmapped Within KeepPairs --outMultimapperOrder Random --outSAMmultNmax 1 --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${sampleID}/${sampleID}.randomonemap. --genomeDir ${genomeIndex} --readFilesIn ${sampleID}/${sampleID}.trimmed_1.fastq.gz ${sampleID}/${sampleID}.trimmed_2.fastq.gz
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// mark duplicate reads in STAR bam
//
//process markduplicates {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus params.markDuplicatesCpus
//
//    tag "Running markduplicates sorted bams using picard"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from bam_ch
//
//    output:
//    file("${sampleID}") into markdupbam1_ch
//    file("${sampleID}") into markdupbam2_ch
//    file("${sampleID}") into markdupbam4_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    java -jar -Xmx40G ${params.picardExecutable} MarkDuplicates I=${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam O=${sampleID}/${sampleID}.markdup.bam M=${sampleID}/${sampleID}.markdup.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT
//
//    java -jar -Xmx40G ${params.picardExecutable} MarkDuplicates I=${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam O=${sampleID}/${sampleID}.rmdup.bam M=${sampleID}/${sampleID}.rmdup.txt REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT
//
//    samtools index ${sampleID}/${sampleID}.markdup.bam
//
//    samtools index ${sampleID}/${sampleID}.rmdup.bam
//
//    samtools sort -n ${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam -o ${sampleID}/${sampleID}.markdup.namesorted.bam
//
//    samtools sort -n ${sampleID}/${sampleID}.randomonemap.Aligned.sortedByCoord.out.bam -o ${sampleID}/${sampleID}.rmdup.namesorted.bam
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// run bamtools stats on marked bams
//
//process bamtools {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//
//    tag "Running bamtools on markdup bams"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from markdupbam1_ch
//
//    output:
//    file("${sampleID}") into bamtools_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    bamtools stats -in ${sampleID}/${sampleID}.markdup.bam > ${sampleID}/${sampleID}.markdup.stats.txt
//
//    sleep ${params.sleepTimeEnd}
//    """
//
//}
//
//// run picard stats
//
//process collect_insert_size_metrics {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Running picard collect insert size metrics on markdup bams"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from bamtools_ch
//
//    output:
//    file("${sampleID}") into picardmetrics_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    java -jar -Xmx40G ${params.picardExecutable} CollectInsertSizeMetrics I=${sampleID}/${sampleID}.markdup.bam H=${sampleID}/${sampleID}.markdup.picardstats.out O=${sampleID}/${sampleID}.markdup.picardstats.txt VALIDATION_STRINGENCY=SILENT
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// make bamtools report
//
//process run_bamtools_report {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Collating results from bamtools"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file dummy1 from picardmetrics_ch.collect()
//
//    output:
//    file("*.txt") into bamtoolsreport_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")
//
//    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir raw_fastqc
//    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir bamtools
//    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir collect_insert_size_metrics 
//
//    ${params.projectDir}/${params.srcDir}/collate_bamtools_and_picard_output.pl ${params.dataDir}/targets.csv ${params.projectDir} ${launchDir}/${params.outputDir} \$linkDir
//
//    rm -r \$linkDir
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// count reads for all genes and transcripts using featureCounts
//
//process featurecounts {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus params.featureCountsJobCpus 
//
//    tag "Running featurecounts for all genes and transcripts"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file sampleID from markdupbam2_ch
//    file gtf from gtf3_ch
//
//    output:
//    file("${sampleID}") into featurecounts_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g gene_id -a ${gtf} -o ${sampleID}/${sampleID}.markdup.genecount.txt ${sampleID}/${sampleID}.markdup.bam
//
//    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g gene_id -a ${gtf} -o ${sampleID}/${sampleID}.rmdup.genecount.txt ${sampleID}/${sampleID}.rmdup.bam
//
//    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g transcript_id -a ${gtf} -o ${sampleID}/${sampleID}.markdup.transcriptcount.txt ${sampleID}/${sampleID}.markdup.bam
//
//    featureCounts -T ${params.featureCountsJobCpus} -O -p -F GTF -t exon -g transcript_id -a ${gtf} -o ${sampleID}/${sampleID}.rmdup.transcriptcount.txt ${sampleID}/${sampleID}.rmdup.bam
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// make count report
//
//process run_featurecounts_report {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Collating results from featureCounts"
//    publishDir path:{params.outputDir},mode: 'symlink'
//
//    input:
//    file dummy from featurecounts_ch.collect()
//    file gtf from gtf4_ch
//
//    output:
//    file("all*.txt") into featurecountsreport_ch
//
//    script:
//    """
//    sleep ${params.sleepTimeStart}
//
//    linkDir=\$(mktemp -d ci-XXXXXXXXXX --tmpdir="${params.projectDir}")
//
//    bash ${params.projectDir}/${params.srcDir}/uber_copy.sh ${launchDir}/trace.txt ${launchDir}/work ${launchDir}/${params.outputDir} \$linkDir featurecounts
//
//    ${params.projectDir}/${params.srcDir}/collate_featurecounts_output.pl ${params.dataDir}/targets.csv ${params.projectDir} ${launchDir}/${params.outputDir} ${launchDir}/${params.resourcesDir} ${gtfName} \$linkDir
//
//    rm -r \$linkDir
//
//    sleep ${params.sleepTimeEnd}
//    """
//}
//
//// make release directory
//
//process make_release {
//    container 'munozcriollojj/nf-pipeline-test:latest'
//    cpus 1
//
//    tag "Making release directory"
//    publishDir path:{params.releaseDir},mode: 'symlink'
//
//    input:
//    file dummy01 from bamtoolsreport_ch.collect()
//    file dummy02 from featurecountsreport_ch.collect()
//    file dummy03 from markdupbam4_ch.collect()
//    file dummy04 from multiqc1_ch
//
//    output:
//    file("*") into release_ch
//
//    script:
//    """
//
//    ${params.projectDir}/${params.srcDir}/make_release.pl ${params.projectDir} ${params.analysisID} targets.csv ${params.dataDir} ${launchDir}/${params.outputDir} ${launchDir}/${params.resourcesDir} ${params.gtfName}
//
//
//    """
//}
