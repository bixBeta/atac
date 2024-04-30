process MACS2 {
        
        maxForks 8
        tag "$id"
        label "process_medium"

        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*narrowPeak"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_peaks.xls"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_summits.bed"

        

        input:

            tuple val(id), path(atac_bam)
            tuple val(id2), path(bg_control)
            val(qval)
            val(fecutoff)
            val(gsize) 

        output:

            tuple val(id), path("*narrowPeak")                  , emit: "narrow_peaks"
            tuple val(id), path("*_peaks.xls")                  , emit: "peaks_xls"
            tuple val(id), path("*_summits.bed")                , emit: "summits_bed"
        
        script:

            """
                macs2 callpeak -t ${atac_bam} \\
                    -f BAMPE \\
                    -n ${id} \\
                    -g ${gsize} \\
                    -q ${qval} \\
                    --nomodel --shift 37 --ext 73 \\
                    --fe-cutoff ${fecutoff} \\
                    --keep-dup all \\
                    --nolambda \\
                    -c ${bg_control}

            """

}


process MACS2ALL {
        
        maxForks 8
        tag "allSamplesMergedPeakset"
        label "process_medium"

        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*narrowPeak"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_peaks.xls"
        publishDir "MACS2_peaks",            mode: "symlink", overwrite: true, pattern: "*_summits.bed"

        

        input:

            val(atac_bam)
            tuple val(id2), path(bg_control)
            val(qval)
            val(fecutoff)
            val(gsize) 

        output:

            path("*narrowPeak")                                 , emit: "all_narrow_peaks"
            path("*_peaks.xls")                                 , emit: "all_peaks_xls"
            path("*_summits.bed")                               , emit: "all_summits_bed"
        
        script:

            b = atac_bam.join(' ')
            """
                macs2 callpeak -t ${b} \\
                    -f BAMPE \\
                    -n allSamplesMergedPeakset \\
                    -g ${gsize} \\
                    -q ${qval} \\
                    --nomodel --shift 37 --ext 73 \\
                    --fe-cutoff ${fecutoff} \\
                    --keep-dup all \\
                    --nolambda \\
                    -c ${bg_control}

            """

}