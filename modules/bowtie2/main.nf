aligner = params.bowtie2
gkey    = params.genome

process BOWTIE2 {

    maxForks 6
    tag "$id"
    label "process_high"

    publishDir "primary_BAMS", mode: "symlink", overwrite: true, pattern: "*.bam*"
    publishDir "STATS",        mode: "symlink", overwrite: true, pattern: "*stat*"
    publishDir "STATS",        mode: "symlink", overwrite: true, pattern: "*log"

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

process MTBLKDUP {

    maxForks 8
    tag "$id"
    label "process_medium"

    publishDir "STATS",             mode: "symlink", overwrite: true, pattern: "*stat*"
    publishDir "DEDUP_BAMS",        mode: "symlink", overwrite: true, pattern: "*DEDUP.bam*"

    input:
        tuple val(id), path(primary_bam)
        tuple val(id), path(primary_bai)
        path(blacklist)

    output:
        tuple val(id), path("*.noMT.bam")                           , emit: "nomt_bam"
        tuple val(id), path("*.noMT.bam.bai")                       , emit: "nomt_bai"
        path("*.noMT.flagstat")                                     , emit: "nomt_flagstat"
        path("*.noMT.idxstats")                                     , emit: "nomt_idxstats"
        
        tuple val(id), path("*.noMT.noBL.bam")                      , emit: "nomt_nobl_bam"
        tuple val(id), path("*.noMT.noBL.bam.bai")                  , emit: "nomt_nobl_bai"
        path("*.noMT.noBL.flagstat")                                , emit: "nomt_nobl_flagstat"
        path("*.noMT.noBL.idxstats")                                , emit: "nomt_nobl_idxstats"

        tuple val(id), path("*.noMT.noBL.dupMarked.bam")            , emit: "dupmarked_bam"
        tuple val(id), path("*.noMT.noBL.dupMarked.bam.bai")        , emit: "dupmarked_bai"
        path("*.noMT.noBL.dupMarked.flagstat")                      , emit: "nomt_nobl_dupmarked_flagstat"
        path("*.noMT.noBL.dupMarked.idxstats")                      , emit: "nomt_nobl_dupmarked_idxstats"        

        path("*.MarkDuplicates.metrics.txt")                        , emit: "dup_stats"

        tuple val(id), path("*DEDUP.bam")                            , emit: "dedup_bam"
        tuple val(id), path("*DEDUP.bam.bai")                        , emit: "dedup_bai"
        path("*.DEDUP.flagstat")                                     , emit: "dedup_flagstat"
        path("*.DEDUP.idxstats")                                     , emit: "dedup_idxstats"

    script:

        """
            samtools view -H ${primary_bam} | cut -f2 | grep "SN:" |  cut -d ":" -f2 | grep -v "MT\\|_\\|\\." | xargs samtools view -b ${primary_bam} > ${id}.noMT.bam

            samtools index ${id}.noMT.bam
            samtools flagstat ${id}.noMT.bam > ${id}.noMT.flagstat
            samtools idxstats ${id}.noMT.bam > ${id}.noMT.idxstats


            bedtools intersect -v -a ${id}.noMT.bam -b ${blacklist} > ${id}.noMT.noBL.bam

            samtools index ${id}.noMT.noBL.bam
            samtools flagstat ${id}.noMT.noBL.bam > ${id}.noMT.noBL.flagstat
            samtools idxstats ${id}.noMT.noBL.bam > ${id}.noMT.noBL.idxstats

            java -jar /myBin/picard.jar \\
                    MarkDuplicates \\
                    INPUT=${id}.noMT.noBL.bam \\
                    OUTPUT=${id}.noMT.noBL.dupMarked.bam \\
                    ASSUME_SORTED=true \\
                    REMOVE_DUPLICATES=false \\
                    METRICS_FILE=${id}.MarkDuplicates.metrics.txt \\
                    VALIDATION_STRINGENCY=LENIENT \\
                    TMP_DIR=tmp

            
            sed -i 's/${id}.noMT.noBL.bam/${id}.noMT.noBL.dupMarked.bam/g' ${id}.MarkDuplicates.metrics.txt   


            samtools index ${id}.noMT.noBL.dupMarked.bam
            samtools flagstat ${id}.noMT.noBL.dupMarked.bam > ${id}.noMT.noBL.dupMarked.flagstat
            samtools idxstats ${id}.noMT.noBL.dupMarked.bam > ${id}.noMT.noBL.dupMarked.idxstats


            samtools view -b -h -F 0X400 ${id}.noMT.noBL.dupMarked.bam > ${id}.DEDUP.bam

            samtools index ${id}.DEDUP.bam
            samtools flagstat ${id}.DEDUP.bam > ${id}.DEDUP.flagstat
            samtools idxstats ${id}.DEDUP.bam > ${id}.DEDUP.idxstats

        """

}

