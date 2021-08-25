# Scan_coverage_genes
Calculating mean coverage &ge; 10x coverage depth for each coding gene of a genome

 <h3>1. Preparing a pileup file for Scan_coverage_genes</h3>
 <code> bedtools bamtobed -i file.sorted.bam > file.sorted.bed</code>
 Then compute depth coverage at each position of the genome:
 <code> bedtools genomecov -i file.sorted.bed -g xxx.txt -d > file.sorted_couv.txt</code>
 
 in progress...
