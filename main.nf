nextflow.enable.dsl=2

// Project Params:
params.sheet            = "sample-sheet.csv"
params.reads            = "$workDir/fastqs/*_*{1,2}.f*.gz"

// Module Params:
params.help             = false
params.listGenomes      = false
params.bwa              = false
params.bowtie2          = false
params.fastp            = false

// Default Params:
params.mode             = "PE"
params.id               = "TREx_ID"
params.fecutoff         = 5
params.qval             = 0.05
params.genome           = null
params.bg               = ''


// Command Line Channels     ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~ 

ch_pin          =    channel.value(params.id)
ch_sheet        =    channel.fromPath(params.sheet)
ch_qval         =    channel.value(params.qval)
ch_fe           =    channel.value(params.fecutoff)
ch_mqc_conf     =    channel.fromPath("${projectDir}/multiqc_config.yaml")
ch_mqc_logo     =    channel.fromPath("${projectDir}/img/trex-mini.png")



// ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~  ~ ~ ~ ~ ~ ~ 

if( params.help ) {

log.info """
A  T  A  C  -  S  E  Q      W  O  R  K  F  L  O  W  -  @bixBeta
=======================================================================================================================================================================
Usage:
    nextflow run https://github.com/bixbeta/atac -r main < args ... >

Args:
    * --help           : Prints this help documentation
    * --listGenomes    : Get extended list of genomes available for this pipeline
    * --id             : TREx Project ID 
    * --sheet          : sample-sheet.csv < default: looks for a file named sample-sheet.csv in the project dir >

        -------------------------------------------
        Sample Sheet Example:    
        label   fastq1          fastq2
        SS1     SS1_R1.fastq.gz SS1_R2.fastq.gz
        SS2     SS2_R1.fastq.gz SS2_R2.fastq.gz  
        .
        .
        . etc.
        -------------------------------------------
    * --bwa         : deprecated
    * --bowtie2     : invokes bowtie2 alignment step
    * --fastp       : invokes fastp trimming
    * --mode        : PE (default)
    * --fecutoff    : 5 (default)
    * --qval        : 0.05 (default)
    * --genome      : reference genome ( Available genomes: hg38, mm10 & dm6 : Use --listGenomes for more details )
    * --bg          : Input Control Sample ID 

"""

    exit 0
}

// BT2 Indices MAP

genomeDir = [

mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/bowtie2index/mm10",
hg38            :"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/BT2.ENSEMBL_INDEX/",
dm6             :"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Bowtie2.Index/dm6"

]

gAlias = [

mm10            :"mouse",
hg38            :"human",
dm6             :"fly"

]

gtfs = [
    
mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/Mus_musculus.GRCm38.96.gtf",
hg38            :"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/Homo_sapiens.GRCh38.96.gtf",
dm6             :"/workdir/genomes/Drosophila_melanogaster/dm6/ENSEMBL/Drosophila_melanogaster.BDGP6.32.106.gtf"

]

blkList = [

mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/no-chr-mm10-blacklist.v2.bed",
hg38            :"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/no-chr-hg38-blacklist.v2.bed",
dm6             :"/workdir/tools/blacklists/dm6-blacklist.v2.bed"

]

gSize = [

mm10            :"mm",
hg38            :"hs",
dm6             :"dm"

]

bwgSize = [

mm10            :2652783500,
hg38            :2913022398,
dm6             :142573017

]


if( params.listGenomes) {
    
    println("")
    log.info """
    Available BT2 Indices
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    genomeDir.each(printMap)

    log.info """
    gtfs
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    gtfs.each(printMap)

    log.info """
    blackLists
    =========================================================================================================================
    """
    .stripIndent()

    printMap = { a, b -> println "$a ----------- $b" }
    blkList.each(printMap)

    exit 0
}



// Load all modules 

include {    FASTP     } from './modules/fastp'
include {    BOWTIE2   } from './modules/bowtie2'
include {    MTBLKDUP  } from './modules/bowtie2'
include {    TAGDIR    } from './modules/homer'
include {    MACS2     } from './modules/macs2'
include {    MACS2ALL  } from './modules/macs2'
include {    FRIP      } from './modules/macs2'
include {    BIGWIG    } from './modules/deeptools'
include {    RSYNC     } from './modules/deeptools'
include {    MQC       } from './modules/multiqc'
include {    MQC2      } from './modules/multiqc'
 



// BOWTIE2 WORKFLOW: 

workflow BTPAIRED {
  
    meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("fastqs/" + row.fastq1), file("fastqs/" + row.fastq2)]] }
                |  view

    
    if( params.fastp ){

        FASTP(meta_ch)
            .set { fastp_out }
    
    } else {

        fastp_out = meta_ch
    
    }

    genome = genomeDir[params.genome]
    ch_genome   = channel.value(genome)

    ch_alias    = channel.value(gAlias[params.genome])
    ch_gtf      = channel.value(gtfs[params.genome])
    ch_blkList  = channel.value(blkList[params.genome])
    ch_gsize    = channel.value(gSize[params.genome])
    ch_bwg_size = channel.value(bwgSize[params.genome])

    if( params.bowtie2){

        BOWTIE2(fastp_out, ch_genome)

        MTBLKDUP(BOWTIE2.out.primary_sorted_bam, 
                 BOWTIE2.out.primary_sorted_bai,
                 ch_blkList)
        
        TAGDIR(MTBLKDUP.out.dedup_bam)

        BIGWIG(MTBLKDUP.out.dedup_bam, MTBLKDUP.out.dedup_bai, ch_bwg_size)

        // RSYNC(BIGWIG.out.big_wig)

        // RSYNC.out.ftp_path.collect().view()

        
        ch_dedup_bams = MTBLKDUP.out.dedup_bam
                         .view()
    
        if(params.bg){

            ch_bg = ch_dedup_bams
                    .filter{ it[0] == params.bg}
                    .collect()
                    .view()

        } else {

            ch_bg = ['', '/']
        }


        //ch_bg_val = channel.value(ch_bg)

        ch_atac_bams = ch_dedup_bams
                            .filter{ it[0] != params.bg}
                            .view()


        MACS2(ch_atac_bams, ch_bg , ch_qval, ch_fe, ch_gsize)


        ch_atac_bams_all = ch_atac_bams
                                .map({ it -> it[1]}).collect()
                                .view()



        MACS2ALL(ch_atac_bams_all, ch_bg , ch_qval, ch_fe, ch_gsize, ch_genome, ch_gtf)

        ch_saf = MACS2ALL.out.saf.collect().view()

        FRIP(ch_atac_bams, ch_saf)

        ch2_mqc = FRIP.out.raw_counts_summary
            .collect()
            .view()

        MQC2(ch2_mqc, ch_mqc_conf, ch_mqc_logo)
    

        


    }

    if( params.genome != null ){

        ch1_mqc = BOWTIE2.out.primary_log
                        .concat(BOWTIE2.out.primary_flagstat)
                        .concat(BOWTIE2.out.primary_idxstats)
                        // .concat(MTBLKDUP.out.nomt_flagstat)
                        // .concat(MTBLKDUP.out.nomt_idxstats)
                        .concat(MTBLKDUP.out.nomt_nobl_dupmarked_flagstat)
                        .concat(MTBLKDUP.out.nomt_nobl_dupmarked_idxstats)
                        .concat(MTBLKDUP.out.dup_stats)
                        .concat(MTBLKDUP.out.dedup_flagstat)
                        .concat(MTBLKDUP.out.dedup_idxstats)
                        .collect()
                        .view()


        MQC(ch1_mqc, ch_mqc_conf, ch_mqc_logo)


    }

}


workflow {

    if (params.mode == "PE" ) {

        BTPAIRED()

    } else {


        error "Invalid Workflow"
        exit 0 
    }


}