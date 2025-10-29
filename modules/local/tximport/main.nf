/*
 * This process takes the kallisto quant output as input and generates the gene expressio nmatrix.
 */
process TXIMPORT {
    tag "tximport_all_samples"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://pdichiaro/lmdseq:latest'

    input:
    path all_folders   // this will be a list of kallisto result folders
    path gtf                            
    path reference                       

    output:
    path("EX_reads_RAW.txt"), emit: matrix   // Gene expression matrix
    path("versions.yml"), emit: versions

    script:
    """
    # Create space-separated list of directories for R script
    DIRS="${all_folders.join(' ')}"

    # Run the R script to generate the expression matrix
    Rscript ${projectDir}/bin/create_reference_db.R \\
        -g ${gtf} \\
        -r ${reference} \\
        -i "\$DIRS" \\
        -o EX_reads_RAW.txt

    # Save versions
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        Rscript: \$(Rscript --version | sed -e "s/Rscript (R) //g")
    END_VERSIONS
    """
}
