/*
 * DESeq2 normalization and quality control analysis
 */
process DESEQ2_NORM_QC {
    tag "deseq2_norm_qc"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker://pdichiaro/lmdseq:latest'

    input:
    path gene_matrix        // Gene expression matrix from TXIMPORT
    path sample_metadata    // Sample metadata file (Sample.txt)

    output:
    path("6_Norm_folder/**"),           emit: norm_files
    path("7_Counts_folder/**"),         emit: count_files  
    path("8_Quality_folder/**"),        emit: qc_files
    path("versions.yml"),               emit: versions

    script:
    """
    # Run DESeq2 normalization and QC
    Rscript ${projectDir}/bin/deseq2_norm_qc.R \\
        --input ${gene_matrix} \\
        --metadata ${sample_metadata} \\
        --outdir . \\
        --min_reads ${params.min_reads}

    # Save versions
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        Rscript: \$(Rscript --version | sed -e "s/Rscript (R) //g")
        DESeq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p 6_Norm_folder/Read_Distribution
    mkdir -p 7_Counts_folder
    mkdir -p 8_Quality_folder
    
    touch 6_Norm_folder/Read_Distribution/Read_Distribution_Raw.pdf
    touch 6_Norm_folder/Read_Distribution/Read_Distribution_Norm_Filt.pdf
    touch 7_Counts_folder/EX_reads_RAW_filt.txt
    touch 7_Counts_folder/EX_reads_NORM_filt.txt
    touch 7_Counts_folder/rLog_reads.txt
    touch 7_Counts_folder/Normalisation_Parameters.txt
    touch 8_Quality_folder/Heatmap_sampleTosample_distances_vstTransformed.pdf
    touch 8_Quality_folder/Heatmap_Poisson_Distance.pdf
    touch 8_Quality_folder/Heatmap_sampleTosample_correlation.pdf
    touch 8_Quality_folder/PCA_rlogTransformedID.pdf
    touch 8_Quality_folder/Percent_var_PCA.pdf
    touch 8_Quality_folder/All_pca.pdf

    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        Rscript: \$(echo "4.3.0")
        DESeq2: \$(echo "1.42.0")
    END_VERSIONS
    """
}