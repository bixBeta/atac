mqcgenome =  params.genome 

process MQC {

    label 'process_mqc'
    
    publishDir "Reports", mode: "symlink", overwrite: true
    
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
    
    publishDir "Reports", mode: "symlink", overwrite: true

    input:

        path "*"
        path(conf)
        path(logo)             

    output:
        path "*html"                    , emit: mqc_out2 

    when:
        
    script:

    """
       export  MQC_GENOME=${mqcgenome} 
       multiqc -n ${params.id}_FRIP.multiqc.report --config ${conf} \\
            --cl-config "custom_logo: ${logo}" -m featureCounts \\
            -b "Please note that the featureCounts M Assigned Column refers to Fragments and Not Reads" .

    """

}