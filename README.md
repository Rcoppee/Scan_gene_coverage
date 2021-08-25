# Scan_coverage_genes
Calculating mean coverage and &ge; 10x coverage depth for each coding gene of a genome. DESCRIPTION IN PROGRESS
We supposed you produced a file_sorted.bam file (for example using Samtools).
<br><br>
 <h3>1. Preparing a pileup file for Scan_coverage_genes analysis</h3>
 <p>Prior to compute depth coverage at each position of the genome, you must convert your sorted.bam file into a bed file with bedtools bamtobed:</p>
 <p><code> bedtools bamtobed -i file_sorted.bam > file_sorted.bed</code></p>
 <p>To compute depth coverage at each position of the genome with bedtools genomecov, you must specified a text file containing the list of chromosomes and corresponding length concomitantly with the bed file. A list_chromosomes.txt file for the Plasmodium falciparum species (version 39 on PlasmoDB) is provided in the data directory.</p>
 <p><code> bedtools genomecov -d -i file_sorted.bed -g list_chromosomes.txt  > file_coverage.txt</code></p>
 <br>
 <h3>2. Calculating mean coverage of each coding gene and percentage of coding gene covered at &ge; 10x depth</h3>
 <p>in progress...</p>
 <br>
 <h3>3. Citation</h3>
 <p>iIf you use this program for your own work, please cite:</p>
 <code>xxxxxxx</code>
