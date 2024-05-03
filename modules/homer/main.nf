process TAGDIR {

    maxForks 8 
    tag "$id"
    label "process_medium"

    publishDir "TAG_DIR",            mode: "symlink", overwrite: true, pattern: "*.tag.dir"


    input:

        tuple val(id), path(dedup_bam)


    output:

        tuple val(id), path("*.tag.dir")                 , emit: "tag_dir"
        path("${id}.tag.dir/*.ucsc.bg.gz")                              , emit: "bedgraph"

    script:

        """
            makeTagDirectory ${id}.tag.dir ${dedup_bam}
            makeUCSCfile ${id}.tag.dir -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq
            
            cd ${id}.tag.dir
            zcat *.ucsc.bedGraph.gz | awk '{if(NR>1) print "chr"\$0; else print \$0}' | gzip > `basename *.ucsc.bedGraph.gz .ucsc.bedGraph.gz`.ucsc.bg.gz
            cd ..
        
        """

}

