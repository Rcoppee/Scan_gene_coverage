# Scan_coverage_genes
Calculating mean coverage &ge; 10x coverage depth for each coding gene of a genome

 <h3>1. Preparing a pileup file for Scan_coverage_genes</h3>
 <p>Prior to compute depth coverage at each position of the genome, you must convert your sorted.bam file into a bed file with bedtools bamtobed:</p>
 <p><code> bedtools bamtobed -i file.sorted.bam > file.sorted.bed</code></p>
 <p>To compute depth coverage at each position of the genome with genomecov, you must specified a text file containing the list of chromosomes and corresponding length concomitantly with the bed file:</p>
 <p><code> bedtools genomecov -i file.sorted.bed -g xxx.txt -d > file.sorted_couv.txt</code></p>
 <p>For example, the list_chromosomes.txt file for the Plasmodium falciparum species may be like that:</p>
 <p><code>CHR1 LENGTH
 
 
 CHR2 LENGTH</code></p>
 in progress...
