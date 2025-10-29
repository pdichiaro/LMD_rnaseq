#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { MULTIQC    } from '../modules/nf-core/multiqc/main'


//
// SUBWORKFLOW: Consisting of local/modules
//
include { FASTQ_FASTQC_TRIMGALORE } from '../subworkflows/local/fastq_fastqc_trimgalore/main' 
include { KALLISTO } from '../subworkflows/local/kallisto/main' 


//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
include { paramsSummaryMap                 } from 'plugin/nf-schema'
include { samplesheetToList                } from 'plugin/nf-schema'
include { paramsSummaryMultiqc             } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML           } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { checkSamplesAfterGrouping        } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { multiqcTsvFromList               } from '../subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness'
include { biotypeInGtf                     } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { getInferexperimentStrandedness   } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'
include { methodsDescriptionText           } from '../subworkflows/local/utils_nfcore_rnaseq_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Header files for MultiQC
//ch_pca_header_multiqc           = file("$projectDir/workflows/assets/multiqc/deseq2_pca_header.txt", checkIfExists: true)
//sample_status_header_multiqc    = file("$projectDir/workflows/assets/multiqc/sample_status_header.txt", checkIfExists: true)
//ch_clustering_header_multiqc    = file("$projectDir/workflows/assets/multiqc/deseq2_clustering_header.txt", checkIfExists: true)
//ch_biotypes_header_multiqc      = file("$projectDir/workflows/assets/multiqc/biotypes_header.txt", checkIfExists: true)
//ch_dummy_file                   = ch_pca_header_multiqc
ch_multiqc_config = Channel.empty()
ch_multiqc_custom_config = Channel.empty()
ch_multiqc_logo = Channel.empty()
ch_replace_names = Channel.empty()
ch_sample_names = Channel.empty()

workflow RNASEQ {
    
    take:
    //ch_samplesheet        // channel: path(sample_sheet.txt)
    ch_versions             // channel: [ path(versions.yml) ]
    ch_fasta                // channel: path(genome.fasta)
    ch_gtf                  // channel: path(gtf)
    ch_transcript_fasta     // channel: path(transcript.fasta)
    ch_index                // channel: [ meta, path(kallisto/index/) ]

    main:
    ch_multiqc_files = Channel.empty()
    ch_versions = Channel.empty()

  

    // -----------------------
    // Run FASTQ preprocessing
    // -----------------------

    // Create channel from input file provided through ch_samplesheet
    ch_samplesheet = Channel.fromPath(params.input)
    ch_input = ch_samplesheet
        .splitCsv(header:true, sep:'\t')
        .map { row ->
            def fastq_dir = file(row.FASTQ_FOLDER)
            def fastq_files = fastq_dir.listFiles().findAll { it.name.endsWith('.fastq.gz') }
            def is_single_end = params.seq_mode == 'SE'

            def r1_files = fastq_files.findAll { it.name.contains('_R1_') }
            def r2_files = fastq_files.findAll { it.name.contains('_R2_') }

            def meta = [
                //id: row.SID,
                id: "${row.SHORT_NAME}_${row.REPLICATE}",
                short_name: row.SHORT_NAME,
                replicate: row.REPLICATE,
                single_end: is_single_end,
                strandedness: row.strandedness
            ]

            def files = is_single_end ? r1_files : r1_files + r2_files
            [ meta, files ]
        }


    // -----------------------
    // Run Subworkflow: cat FASTQs, rename files, read QC and trim adapters
    // -----------------------

    if (params.trimmer == 'trimgalore') {
        FASTQ_FASTQC_TRIMGALORE (
            ch_input,
            params.skip_fastqc,
            params.skip_trimming
        )
    }
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_FASTQC_TRIMGALORE.out.multiqc_files)
    ch_versions = ch_versions.mix(FASTQ_FASTQC_TRIMGALORE.out.versions)
    ch_trimmed_reads = FASTQ_FASTQC_TRIMGALORE.out.reads
    
    //ch_trimmed_reads.collect { println it } // This will give you the current contents of the channel


    // -----------------------
    // Run Subworkflow: Alignment with KALLISTO, generate pseudo-bam, raw count matrix, and DESeq2 normalization/QC
    // -----------------------
    ch_chromosomes = params.chromosomes ? Channel.fromPath(params.chromosomes) : Channel.empty()        
    if (params.aligner == 'kallisto') {
        KALLISTO(
            ch_trimmed_reads,
            ch_index,
            ch_gtf,
            ch_chromosomes,
            params.kallisto_quant_fraglen,
            params.kallisto_quant_fraglen_sd,
            params.bin_size,
            params.reference ? file(params.reference) : [],
            ch_samplesheet
        )
    }
    ch_multiqc_files = ch_multiqc_files.mix(KALLISTO.out.multiqc_files)
    ch_versions = ch_versions.mix(KALLISTO.out.versions)
    ch_bigwig = KALLISTO.out.bigwig
    ch_raw_matrix = KALLISTO.out.matrix
    ch_qc_files = KALLISTO.out.qc_files


    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'LMD_rnaseq_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_report = Channel.empty()

    if (!params.skip_multiqc) {

        // Load MultiQC configuration files
        ch_multiqc_config        = Channel.fromPath("$projectDir/workflows/assets/multiqc/multiqc_config.yml", checkIfExists: true)
        ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()
        ch_multiqc_logo          = params.multiqc_logo   ? Channel.fromPath(params.multiqc_logo)   : Channel.empty()

        // Prepare the workflow summary
        ch_workflow_summary = Channel.value(
            paramsSummaryMultiqc(
                paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
            )
        ).collectFile(name: 'workflow_summary_mqc.yaml')

        // Prepare the methods section
        ch_methods_description = Channel.value(
            methodsDescriptionText(
                params.multiqc_methods_description
                    ? file(params.multiqc_methods_description)
                    : file("$projectDir/workflows/assets/multiqc/methods_description_template.yml", checkIfExists: true)
            )
        ).collectFile(name: 'methods_description_mqc.yaml')

        // Add summary, versions, and methods to the MultiQC input file list
        ch_multiqc_files = ch_multiqc_files
            .mix(ch_workflow_summary)
            .mix(ch_collated_versions)
            .mix(ch_methods_description)

        // Provide MultiQC with rename patterns to ensure it uses sample names
        // for single-techrep samples not processed by CAT_FASTQ, and trims out
        // _raw or _trimmed

        ch_name_replacements = ch_input
            .map{ meta, reads ->
                def name1 = file(reads[0][0]).simpleName + "\t" + meta.id + '_1'
                def fastqcnames = meta.id + "_raw\t" + meta.id + "\n" + meta.id + "_trimmed\t" + meta.id
                if (reads[0][1] ){
                    def name2 = file(reads[0][1]).simpleName + "\t" + meta.id + '_2'
                    def fastqcnames1 = meta.id + "_raw_1\t" + meta.id + "_1\n" + meta.id + "_trimmed_1\t" + meta.id + "_1"
                    def fastqcnames2 = meta.id + "_raw_2\t" + meta.id + "_2\n" + meta.id + "_trimmed_2\t" + meta.id + "_2"
                    return [ name1, name2, fastqcnames1, fastqcnames2 ]
                } else{
                    return [ name1, fastqcnames ]
                }
            }
            .flatten()
            .collectFile(name: 'name_replacement.txt', newLine: true)

        MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList(),
            ch_name_replacements,
            []
        )
        ch_multiqc_report = MULTIQC.out.report
    }


    emit:
    multiqc_report      = ch_multiqc_report      // channel: /path/to/multiqc_report.html
    raw_matrix          = ch_raw_matrix          // channel: [ path(EX_reads_RAW.txt) ]
    qc_files            = ch_qc_files            // channel: [ path(quality_control_files) ]
    versions            = ch_versions            // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
