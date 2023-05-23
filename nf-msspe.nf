#!/usr/bin/env nextflow
/*
========================================================================================
                         mfurla/bproject
========================================================================================
 MaestSi/nf-msspe analysis pipeline.
 #### Homepage / Documentation
 https://github.com/MaestSi/nf-msspe
----------------------------------------------------------------------------------------
*/

def helpMessage() {
        log.info"""
    Usage:
    nextflow -c nf-msspe.conf run nf-msspe.nf --fasta_file = "/path/to/file.fasta" --primers_file = "/path/to/primers.fasta" -profile docker

    Mandatory argument:
    -profile                                                              Configuration profile to use. Available: docker, singularity
    Other mandatory arguments which may be specified in the nf-msspe.conf file
    --fasta_file = "/path/to/file.fasta"                                  Path to fasta file with genomes of interest
    --primers_file = "/path/to/primers.fasta"                             Path to output fasta file with primers
    --ovlp_window_size = 250                                              Size of windows overlaps
    --search_window_size = 50                                             Search window size at the extremities of windows
    --kmer_size = 13                                                      Primers length
    --num_acc_miss = 0                                                    Number of accepted missed sequences, to reduce costs
    --num_max_it = 1000                                                   Maximum number of iterations
    --scripts_dir = "/path/to/scripts/dir"                                Path to directory containing MSSPE.R script
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Input of fasta file
Channel
	.fromPath(params.fasta_file, checkIfExists:true)
	.into{fasta_msa}

// From fasta to multiple sequence alignment
process msa {
    input:
		tuple val(fasta) from fasta_msa

    output:
    	val(fasta) into msa_msspe
    script:
    if(params.msa)
    """
    	mkdir -p \$(basename ${params.primers_file}
    	mfa_file=\$(echo \$(basename \$(realpath ${params.fasta_file}) | sed \'s/\\.fa.*/.mfa/\'))
    	mafft --auto --thread ${task.cpus} --adjustdirectionaccurately ${fasta} > \$mfa_file 
    """
	else
	"""
		mfa_file=\$(echo \$(basename \$(realpath ${params.fasta_file}) | sed \'s/\\.fa.*/.mfa/\'))
		ln -s \$mfa_file .
    """
}

// From fasta to multiple sequence alignment
process msspe {
    input:
		tuple val(fasta) from msa_msspe

    output:

    script:
    if(params.msspe)
    """
    	mkdir -p \$(basename ${params.primers_file}
    	mfa_file=\$(echo \$(basename ${params.fasta_file}) | sed \'s/\\.fa.*/.mfa/\'))
    	
    	Rscript ${params.scripts_dir}/MSSPE.R ovlp_window_size=${params.ovlp_window_size} 
    	search_window_size=${params.search_window_size} kmer_size=${params.kmer_size} 
    	num_acc_miss=${params.num_acc_miss} num_max_it={params.num_max_it} 
    	mfa_file=\$mfa_file primers_file=${params.primers_file}
    """
	else
	"""
		echo "Skipped."
    """
}

