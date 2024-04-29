aligner = params.bowtie2
gkey    = params.genome

process BOWTIE2M {

    maxForks 6
    tag "$id"
    label "process_high"

    publishDir "primary_BAMS", mode: "symlink", overwrite: true, pattern: "*.primary.bam"

    input:
        tuple val(id), path(trimmed)
        path genome


    output:
        tuple val(id), path("*.primary.bam")        , emit: "primary_bam"
        tuple val(id), path("*.primary.log")        , emit: "primary_log"


    script:

        """
            (bowtie2 \
            --no-unal \
              -x  ${gkey} \
              -1 ${trimmed[0]} -2 ${trimmed[1]} \
              --threads 12 \
              -S - | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${id}.primary.bam)2>${id}.primary.log

        """


}