# Nextflow Pipeline for ATAC-seq runs on GG02 
![](img/trex-mini.png)

<hr>

This pipeline can be directly used on GG02. To envoke the pipeline please make sure you have nextflow available in your path. 
You may use the [ following guide](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=638#c>) to help you configure nextflow for your account. 


To see if this pipeline works on your account simply run the following command on your terminal:

```
nextflow run https://github.com/bixbeta/atac -r main --help
```

If successful, you may see the following output on your terminal console:

![](img/success.png)


<hr>
It is always a good idea to run the pull command before executing the pipeline. This will ensure that the user is using the latest branch version of the pipeline.
Use the following command to ensure the usage of the latest version of the pipeline:

```
nextflow pull https://github.com/bixbeta/atac -r main 
```
<hr>

<details>
  <summary> DAG 
    ![](DAG/DAG.png)
</summary>
</details>
