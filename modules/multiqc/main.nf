mqcgenome =  params.genome 

process MQC {

    label 'process_mqc'
    
    publishDir "Reports", mode: "move", overwrite: true
    input:

        path "*"
        path(conf)
        path(logo)             

    output:
        path "*html"                    , emit: mqc_out  

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome} 
       multiqc -n ${params.id}_bt2_Primary.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" .

    """

}


process MQC2 {

    label 'process_mqc'
    
    publishDir "Reports", mode: "move", overwrite: true
    input:

        path "*"
        path(conf)
        path(logo)             

    output:
        path "*html"                    , emit: mqc2_out  

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome} 
       multiqc -n ${params.id}_bt2_FRIP.multiqc.report --config ${conf} --cl-config "custom_logo: ${logo}" .

    """

}