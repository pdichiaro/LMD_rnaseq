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
    path reference                      // Optional: can be empty list []                       

    output:
    path("EX_reads_RAW.txt"), emit: matrix   // Gene expression matrix
    path("versions.yml"), emit: versions

    script:
    def reference_arg = reference.toString() != "[]" ? "-r ${reference}" : ""
    """
    # Create space-separated list of directories for R script
    DIRS="${all_folders.join(' ')}"

    # Run the R script to generate the expression matrix
    # Reference is optional - only pass if provided
    Rscript ${projectDir}/bin/create_reference_db.R \\
        -g ${gtf} \\
        ${reference_arg} \\
        -i "\$DIRS" \\
        -o EX_reads_RAW.txt

    # Save versions
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        Rscript: \$(Rscript --version | sed -e "s/Rscript (R) //g")
    END_VERSIONS
    """
}
