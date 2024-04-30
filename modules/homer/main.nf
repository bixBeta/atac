process TAGDIR {

    maxForks 8 
    tag $id
    label "process_medium"

    input:

        tuple val(id), path(dedup_bam)


    output:

        tuple val(id), path("*tag.dir")                 , emit: "tag_dir"
        
    script:

        """
            makeTagDirectory ${id}.tag.dir ${dedup_bam}

        """

}