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
        tuple val(id), path("*.primary.sorted.bam")             , emit: "primary_sorted_bam"
        tuple val(id), path("*.primary.sorted.bam.bai")         , emit: "primary_sorted_bai"
        path("*.primary.log")                                   , emit: "primary_log"
        path("*.primary.flagstat")                              , emit: "primary_flagstat"
        path("*.primary.idxstats")                              , emit: "primary_idxstats"

    script:

        """
            (bowtie2 \
            --no-unal \
              -x  ${genome}/hg38 \
              -1 ${trimmed[0]} -2 ${trimmed[1]} \
              --threads 24 \
              -S - | samtools view -@ 24 -b -h -F 0x0100 -O BAM -o ${id}.primary.bam)2>${id}.primary.log


        samtools sort ${id}.primary.bam > ${id}.primary.sorted.bam
        samtools index ${id}.primary.sorted.bam 

        samtools flagstat ${id}.primary.sorted.bam > ${id}.primary.flagstat
        samtools idxstats ${id}.primary.sorted.bam > ${id}.primary.idxstats

        """


}