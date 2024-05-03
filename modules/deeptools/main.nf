process BIGWIG {

    maxForks 8 
    tag "$id"
    label "process_deeptools"

    publishDir "BIG_WIGS",            mode: "symlink", overwrite: true, pattern: "*.bw"


    input:

        tuple val(id), path(dedup_bam)
        tuple val(id), path(dedup_bai)
        val egsize

    output:

        tuple val(id), path("*.bw")                 , emit: "big_wig"
        
    script:

        """
            bamCoverage -b ${dedup_bam} -o ${id}.coverage.bw \\
                -bs 50 \\
                -p 20 \\
                --normalizeUsing CPM \\
                --effectiveGenomeSize ${egsize} 


        """

}

