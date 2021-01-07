# Parent-of-origin DNA methylation in the bumblebee, Bombus terrestris.

This project is a collaboration between the labs of Dr. Eamonn Mallon (University of Leicester, UK), Prof. Tom Wenseleers (KU Leuven, BE) and Dr. Laura Ross (Unviersity of Edinburgh, UK). We aim to identity parent-of-origin DNA methylation in a species of bumblebee with the prediction that it may act as a mechanism of generating parent-of-origin gene expression, i.e. imprinted genes. 

---

## Task 1: Create alternate reference genomes for each parent using SNPs called from WGS

- Quality check, trim and align genome sequencing to ref genome *DONE*
- Deuplicate and add read groups *DONE*
- Realign around indels *DONE*
- Call SNPs and filter *DONE*
- Make alternate reference genome *DOING*
- Make SNP files for informative queen and male SNPs
- Calculate average coverage per SNP and make a bar plot to show the number of filtered SNPs etc.

See [here](./Genome_processing/genome_processing.md) for more details.

---

## Future Tasks:

### Differential DNA methylation between castes
- This is done (align use the standard Bter ref to make sure the number of Cs is consistent): consider using the data aligned to the alternate reference genomes and making individual C files for weighted methylation calculation
- Pull all methods/results/scripts into GitHub

### Determine parent-of-origin DNA methylation
- Align samples to the alternate reference genomes
- Intersect to pull out reads which have informative SNP information and assign to parental allele bam files
- Check the coverage level for the new parental read sets
- Carry out the differential methylation using methylkit pipeline
- Make some informative plots
- If enough sites, annotate with gene ID
    - Carry out GO enrichment analysis

- Consider using DAMEfinder to confirm the results

### Determine parent-of-origin DNA methylation in the honeybee data for comparison
- Gather all data (parental genomes and WGBS of offspring)
- Run through the above pipeline to generate a list of parent-of-origin sites

### Comparative analyses
- If enough genes with parent-of-origin DNA methylation
    - Compare with previous parent-of-origin expression
    - Compare with honeybee parent-of-origin DNA methylation genes (if any) through orthology database

--- 

For further information about this project contact: hollie_marshall@hotmail.co.uk.

