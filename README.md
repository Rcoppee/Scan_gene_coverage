# Scan_coverage_genes
Calculating mean coverage &ge; 10x coverage depth for each coding gene of a genome

 <h3>1. Preparing a pileup file for Scan_coverage_genes</h3>
 <p>Prior to compute depth coverage at each position of the genome, you must convert your sorted.bam file into a bed file with bedtools bamtobed:</p>
 <p><code> bedtools bamtobed -i file_sorted.bam > file_sorted.bed</code></p>
 <p>To compute depth coverage at each position of the genome with bedtools genomecov, you must specified a text file containing the list of chromosomes and corresponding length concomitantly with the bed file. A list_chromosomes.txt file for the Plasmodium falciparum species (version 39 on PlasmoDB) is provided in the data directory.</p>
 <p><code> bedtools genomecov -d -i file_sorted.bed -g list_chromosomes.txt  > file_coverage.txt</code></p>
 <p></p>
 in progress...
