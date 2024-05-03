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


process RSYNC {
 
    tag "$id"
    label "process_low"

    input:

        tuple val(id), path(big_wig)     

    output:

        path("ftp.path")                                , emit: "ftp_path"
        
    shell:
        
        
        '''
            rsync -av !{big_wig} fa286@cbsugg.biohpc.cornell.edu:/workdir/ftp/fa286/!{params.id}
            echo "ftp://cbsuftp.biohpc.cornell.edu/cbsugg/fa286/!{params.id}/!{big_wig}" > ftp.path

        '''

}