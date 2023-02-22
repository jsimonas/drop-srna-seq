#!/usr/bin/env nextflow
/*
========================================================================================
                                    drop-srna-seq
========================================================================================
 jsimonas/drop-srna-seq analysis pipeline.
 #### homepage / documentation
 https://github.com/jsimonas/drop-srna-seq
----------------------------------------------------------------------------------------
*/

// using DSL-2
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run jsimonas/drop-srna-seq --run_dir 'path/to/bcl_folder' --sample_sheet 'path/to/extended_sample_sheet.xlsx' -profile singularity

    Mandatory arguments:
      --run_dir [path/to/folder]      Path to input data (must be surrounded with quotes)
      --run_module [str]              Pipeline module to run. Can be set as "complete", "demux" or "fastq". If the latter selected, sample sheet is not required. Default: "complete".
      --sample_sheet [file]           Full path to extended sample sheet file. Example can be found at solo-in-drops/assets/extended_sample_sheet_template.xlsx
      --sequencer [str]               Sequencer used to generate the data. Default: "nextseq". Can be set as "nextseq" or "miseq".
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated).
                                      Available: conda, docker and singularity

    References:                       If not specified in the configuration file or you wish to overwrite any of the references
      --star_index [path/to/folder]   Path to star index directory (same as --genomeDir parameter in STAR)
      --barcode_whitelist [file]      Path to cell barcode list (a text file containing one barcode sequence per line)
    
    STARsolo arguments:               If not specified, the default parameters will be used
      --bc_read_length [int]          Read length of cell barcode read. Default: equal to sum of BC + UMI
      --solo_multi_mappers            Allow multi-gene read quantification. Can be set as "Uniform", "PropUnique", "EM", "Rescue" or any combination of these options. Default: "Uniform"

    bcl2fastq arguments:              If not specified, the default parameters will be used
      --barcode_mismatches            Allowed missmaches in library index for demultiplexing. Default: 1
      --write_fastq                   Output raw FASTQ files for each library. Default: false
    
    Other options:
      --outdir [file]                 The output directory where the results will be saved
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic
    
    AWSBatch options:
      --awsqueue [str]                The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]               The AWS Region for your AWS Batch job to run on
      --awscli [str]                  Path to the AWS CLI tool  
    """.stripIndent()
}

// import modules
include { demux_wf } from './modules/demux'


// main workflow
workflow {
    
    // show help message
    if (params.help) {
        helpMessage()
        exit 0
    }
    
    // catch run name
    custom_runName = params.name
    if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
        custom_runName = workflow.runName
    }

    // 
    // validate parameters
    //
        
    // check run directory
    if (params.run_dir){
        runDir = file(params.run_dir, checkIfExists: true)
    } else {
        exit 1, "Input directory not found!"
    }
    runName = runDir.getName()

    // check if correct sequencer is provided
    if (!(params.sequencer.equals('nextseq') || params.sequencer.equals('miseq'))){
        exit 1, "Unsupported sequencer provided! Currently nextseq or miseq are supported"
    }

    // check if STAR index is provided
    if( params.star_index ){
        star_index = Channel
            .fromPath(params.star_index)
            .ifEmpty { exit 1, "STAR index not found: ${params.star_index}" }
    }

    // check if barcode whitelist is provided
    if( params.barcode_whitelist ){
        barcode_whitelist = Channel
            .fromPath(params.barcode_whitelist)
            .ifEmpty { exit 1, "The barcode whitelist not found: ${params.barcode_whitelist}" }
    }

    // check if sample sheet is provided
    if (params.sample_sheet && (params.run_module.equals('complete') || params.run_module.equals('demux'))){
        sheet_file = file(params.sample_sheet, checkIfExists: true)
    } else if (params.run_module.equals('fastq')) {
        sheet_file = Channel.empty()
    } else {
        exit 1, "The extended sample sheet is not provided! The template can be found at drop-srna-seq/assets/misc/extended_sample_sheet_template.xlsx"
    }

    // check if run module is provided correctly
    if (!(params.run_module.equals('complete') || params.run_module.equals('demux') || params.run_module.equals('fastq'))){
        exit 1, "Uncorrect pipeline run module was provided! Can be set as 'complete', 'demux' or 'fastq' module."
    }

    // 
    // main workflow 
    //
    
    demux_wf()

}