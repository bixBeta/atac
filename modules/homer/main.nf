process TAGDIR {

    maxForks 8 
    tag "$id"
    label "process_medium"

    publishDir "TAG_DIR",            mode: "symlink", overwrite: true, pattern: "*.tag.dir"


    input:

        tuple val(id), path(dedup_bam)


    output:

        tuple val(id), path("*tag.dir")                 , emit: "tag_dir"
        
    script:

        """
            makeTagDirectory ${id}.tag.dir ${dedup_bam}
            makeUCSCfile ${id}.tag.dir -o auto -fsize 1e10 -res 1 -color 106,42,73 -style chipseq

        """

}

