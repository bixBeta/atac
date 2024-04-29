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

"""

    exit 0
}

// BT2 Indices MAP

genomeDir = [

mm10            :"/workdir/genomes/Mus_musculus/mm10/ENSEMBL/bowtie2index/mm10",
hg38            :"/workdir/genomes/Homo_sapiens/hg38/ENSEMBL/BT2.ENSEMBL_INDEX/hg38*",
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

include {    FASTPM    } from './modules/fastp'
include {    BOWTIE2M  } from './modules/bowtie2'
include {    MQC       } from './modules/multiqc'
 



// BOWTIE2 WORKFLOW: 

workflow BTPAIRED {
  
    meta_ch = ch_sheet
                |  splitCsv( header:true )
                |  map { row -> [row.label, [file("fastqs/" + row.fastq1), file("fastqs/" + row.fastq2)]] }
                |  view

    
    if( params.fastp ){

        FASTPM(meta_ch)
            .set { fastp_out }
    
    } else {

        fastp_out = meta_ch
    
    }

    genome = genomeDir[params.genome]
    ch_genome   = channel.value(genome)

    ch_alias    = channel.value(gAlias[params.genome])
    ch_gtf      = channel.value(gtfs[params.genome])
    ch_blkList  = channel.value(blkList[params.genome])

    if( params.bowtie2){

        BOWTIE2M(fastp_out, ch_genome)

    }

    if( params.genome != null ){

        ch1_mqc = BOWTIE2M.out.primary_log
                        .concat(BOWTIE2M.out.primary_flagstat)
                        .concat(BOWTIE2M.out.primary_idxstats)
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