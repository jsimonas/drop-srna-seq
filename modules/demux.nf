#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// import the multiqc process, while specifying the location for the outputs
//include { multiqc as multiqc_demux } from './multiqc' addParams(output_subfolder: 'demux')

//convert extended to standard sample sheet
process convert_sample_sheet {
    tag "$sheet_file"
    label 'process_low'
    publishDir path: "${params.outdir}/", mode: 'copy'
 
    input:
    path sheet_file

    when:
    params.run_module.equals('complete') || params.run_module.equals('demux') 

    output:
    path "*.csv", emit: sheet
    
    script:
    """
    convert_to_samplesheet.py --file $sheet_file --out "standard_samplesheet.csv"
    """
}

//convert bcl to fastq files
process bcl_to_fastq {
    tag "$runName"
    label 'process_high'
    publishDir path: "${params.outdir}/", pattern: "*/*.fastq.gz", mode: 'copy',
        saveAs: { filename -> 
            if (params.write_fastq) filename
            else null
        }
    publishDir path: "${params.outdir}/", pattern: "*.fastq.gz", mode: 'copy',
        saveAs: { filename -> 
            if (params.write_fastq) "Undetermined/fastqs/$filename"
            else null
        }

    input:
    path(sheet)
    path(runDir)

    when:
    params.run_module.equals('complete') || params.run_module.equals('demux') 

    output:
    path "*/**{R1,R2,R3}_001.fastq.gz", emit: fastqs
    path "*{R1,R2,R3}_001.fastq.gz", emit: und_fastqs
    path "Stats", emit: demux_stats

    script:
    
    """
    bcl2fastq \\
    --runfolder-dir $runDir \\
    --output-dir . \\
    --sample-sheet $sheet \\
    --mask-short-adapter-reads 0 \\
    --minimum-trimmed-read-length 0 \\
    --use-bases-mask 'y*,I*,y*,y*' \\
    --no-lane-splitting \\
    --create-fastq-for-index-reads \\
    --barcode-mismatches $params.barcode_mismatches \\
    --processing-threads $task.cpus
    """
}

workflow demux_wf{

    take:
    sheet_file
    runDir

    main:

    //convert extended to standard sample sheet
    convert_sample_sheet(sheet_file)

    //convert bcl to fastq files
    bcl_to_fastq(convert_sample_sheet.out.sheet, runDir)

    // combine the stats reports
    //multiqc_demux(bcl_to_fastq.out.demux_stats.toSortedList())
    
    // add project name
    fqname_fqfile_ch = bcl_to_fastq.out.fastqs
    .flatten()
    .map{
        file -> [file.getParent().getName(), file]
        }
    undetermined_fqfile_ch = bcl_to_fastq.out.und_fastqs
    .flatten()
    .map{
        file -> ["Undetermined", file]
    }
    // combine files
    fastqs = Channel.empty()
    fastqs_ch = fastqs.mix(fqname_fqfile_ch, undetermined_fqfile_ch)
    fastqs_ch
    .view()
    
    emit:
    fastq_files = fastqs_ch

}

