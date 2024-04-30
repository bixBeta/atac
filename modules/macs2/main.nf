process MACS2 {
        
        maxForks 8
        tag "$id"
        label "process_medium"

        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*narrowPeak"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_peaks.xls"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_summits.bed"

        

        input:

            tuple val(id), path(atac_bams)
            val(bg_control) from params.input ? Channel.fromPath(params.input) : Channel.empty()

        output:

            tuple val(id), path("*narrowPeak")
            tuple val(id), path("*_peaks.xls")
            tuple val(id), path("*_summits.bed")
        
        script:

}