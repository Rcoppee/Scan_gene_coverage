# Scan_gene_coverage
<p>An algorithm for calculating the coverage depth for each coding gene and the percentage of each gene covered at &ge; 10X depth.<br>
 The ouput file will give a list of gene identifiers with the corresponding mean coverage and their proportion covered at &ge; 10X depth.<br>
 We suppose you produced a <i>sorted.bam</i> file at this step (for example using Samtools).<br>
 The algorithm was written for <i><b>Plasmodium falciparum</b></i> species, but may also be used for other <i><b>Plasmodium</b></i> species and <b>other organisms</b> (by supplying <i>fasta</i> and <i>gff</i> files).</p>
<br>
 <h3>1. Preparing a per-base coverage file for Scan_gene_coverage analysis</h3>
 <p>Prior to compute depth coverage at each position of the genome, you must convert your <i>sorted.bam</i> file into a <i>sorted.bed</i> file with bedtools bamtobed:</p>
 <p><code> bedtools bamtobed -i file_sorted.bam > file_sorted.bed</code></p>
 <p>To compute depth coverage at each position of the genome with bedtools genomecov, you must specified a text file containing the list of chromosomes and corresponding length concomitantly with the <i>sorted.bed</i> file. A <i>list_chromosomes.txt</i> file for the <i>Plasmodium falciparum</i> species (version 39 on PlasmoDB) is provided in the <code>data</code> directory.</p>
 <p><code> bedtools genomecov -d -i file_sorted.bed -g list_chromosomes.txt > file_coverage.txt</code></p>
 <br>
 <h3>2. Calculating mean coverage of each coding gene and percentage of coding gene covered at &ge; 10x depth</h3>
 <p>To calculate the mean coverage of each coding gene and percentage of coding gene covered at &ge; 10x depth, you must provide a <i>species.gff</i> file (that contains the coordinates of each coding gene), a reference genome in <i>fasta</i> format (it must be the same version as the provided <i>gff</i> file), and your <i>file_coverage.txt</i> file previously obtained with bedtools genomecov.</p>
 <p>An example of reference genome in <i>fasta</i> format and corresponding <i>gff</i> file are provided in the <code>data</code> directory.</p>
 <p><code>python3 Scan_gene_coverage.py -p file_coverage.txt -f reference_genome.fasta -g reference_coordinates.gff -o output.txt</code></p>
 <br>
 <h3>3. Citation</h3>
 <p>If you use this program for your own work, please cite:</p>
 <p>Coppée et al. 5WBF: A low-cost and straightforward whole blood filtration method suitable for whole-genome sequencing of <i>Plasmodium falciparum</i> clinical isolates. (2021) In preparation.</i></p>
