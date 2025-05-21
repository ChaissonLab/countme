This is a simple program to count the average methylation of CpG
islands in intervals of reads that map to intervals in a reference.
This way, if there are expansions in the read that are contained
at the interval on the reference, the CpGs in the expansion are counted.

The input bam must be phased with HP tags.
The result that is output is the mean and variance of CpG, as well as the
number of CpGs for each haplotype over a region.

The output columns are:
1-3 chrom start end (from bed flie)
4-7 hap1: nReads avgReadLen avgCpGCount avgMeth
8-11 hap2: nReads avgReadLen avgCpGCount avgMeth

The nreads are the number of reads from the corresponding haplotype that
overlap the interval, and avgReadLen is the length of each read between
the bases that are mapped to the starting and ending positions of the
interval.  This should roughly equal the starting and ending positions of the
interval.  The avgCpG count is according to the read itself rather
than the genome, and the avgMeth is the average of the methylation scores
for all CpG islands for all reads from a haplotype.






