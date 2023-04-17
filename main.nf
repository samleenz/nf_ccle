// a nextflow pipeline

// Written by Sam Lee
// Davis Lab

include { RUN_FASTQC } from "${params.module_dir}/fastqc"
include { RUN_MULTIQC } from "${params.module_dir}/multiqc"
include { RUN_TRIM_READS } from "${params.module_dir}/cutadapt"
include { RUN_STAR_GENOME_GENERATE } from "${params.module_dir}/star-index" 
include { RUN_STAR_ALIGN as RUN_STAR_PASS1 } from "${params.module_dir}/star-align" 
include { RUN_STAR_ALIGN as RUN_STAR_PASS2 } from "${params.module_dir}/star-align" 
include { RUN_FEATURECOUNTS } from "${params.module_dir}/featureCounts" 
include { RUN_FEATURECOUNTS_COMBINE } from "${params.module_dir}/featureCounts-combine" 

workflow {

    // get the fastq file channel from the samplesheet info
    read_pairs_ch = Channel.fromPath(params.samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map(row -> tuple(row.sample, [file(row.read1), file(row.read2)] ))


    // Trim reads with cutadapt
    trimmed_reads_ch = RUN_TRIM_READS(read_pairs_ch)





    // STAR ALIGNMENT -------

    // make star genome index
    star_idx_ch = RUN_STAR_GENOME_GENERATE(
        params.genome, // genome
        params.gtf, // gtf
        params.tx_ver, // ensembl version
        params.read_len // read_length
    )
  

    // pass one
    (star_align_pass1_ch, star_junctions_pass1_ch, star_logs_pass1_ch) = RUN_STAR_PASS1(
        trimmed_reads_ch,
        star_idx_ch,
        params.read_len,
        file(params.juncs)
    )

    // pass two
    (star_align_pass2_ch, star_junctions_pass2_ch, star_logs_pass2_ch) = RUN_STAR_PASS2(
        trimmed_reads_ch,
        star_idx_ch,
        params.read_len,
        star_junctions_pass1_ch.collect()
    )

    // FEATURECOUNTS -------- 

    (fc_counts_ch, fc_logs_ch) = RUN_FEATURECOUNTS(
        star_align_pass2_ch,
        params.gtf
    )

    fc_combine_ch = RUN_FEATURECOUNTS_COMBINE(fc_counts_ch.collect())



    // QC --------------

    fastqc_ch = RUN_FASTQC(trimmed_reads_ch)

    multiqc_ch = RUN_MULTIQC(
        fastqc_ch
            .concat(star_logs_pass2_ch)
            .concat(fc_logs_ch)
            .collect()
        )



}

workflow.onComplete {
    log.info( 
        workflow.success ? 
        "Pipeline Complete!\n" : 
        "Pipeline Failed.\n" )
}