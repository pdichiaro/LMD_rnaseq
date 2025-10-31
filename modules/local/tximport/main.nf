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
    def reference_arg = (reference.toString() != "[]" && reference.toString() != "") ? "-r ${reference}" : ""
    """
    # Validate inputs
    echo "Input validation:"
    echo "GTF file: ${gtf}"
    echo "Reference: ${reference}"
    echo "Number of kallisto folders: ${all_folders.size()}"
    
    # Check GTF file exists
    if [ ! -f "${gtf}" ]; then
        echo "Error: GTF file not found: ${gtf}"
        exit 1
    fi
    
    # Check all kallisto folders exist and contain abundance.h5
    for dir in ${all_folders.join(' ')}; do
        if [ ! -d "\$dir" ]; then
            echo "Error: Kallisto output directory not found: \$dir"
            exit 1
        fi
        if [ ! -f "\$dir/abundance.h5" ]; then
            echo "Error: abundance.h5 not found in \$dir"
            exit 1
        fi
        echo "✓ Valid kallisto output: \$dir"
    done

    # Create space-separated list of directories for R script
    DIRS="${all_folders.join(' ')}"

    # Run the R script to generate the expression matrix
    echo "Running tximport with the following parameters:"
    echo "  GTF: ${gtf}"
    echo "  Reference: ${reference_arg}"
    echo "  Input dirs: \$DIRS"
    echo "  Output: EX_reads_RAW.txt"
    
    Rscript ${projectDir}/bin/create_reference_db.R \\
        -g ${gtf} \\
        ${reference_arg} \\
        -i "\$DIRS" \\
        -o EX_reads_RAW.txt

    # Validate output was created
    if [ ! -f "EX_reads_RAW.txt" ]; then
        echo "Error: Output file EX_reads_RAW.txt was not created"
        exit 1
    fi
    
    # Show output summary
    echo "Output file created successfully:"
    echo "  File: EX_reads_RAW.txt"
    echo "  Size: \$(wc -l < EX_reads_RAW.txt) lines"
    echo "  Columns: \$(head -1 EX_reads_RAW.txt | tr '\t' '\n' | wc -l) columns"

    # Save versions
    cat <<-END_VERSIONS > versions.yml
    "\${task.process}":
        Rscript: \$(Rscript --version 2>&1 | sed 's/R scripting front-end version //')
    END_VERSIONS
    """
}
