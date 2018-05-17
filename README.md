# NGS-JAK2
Investigating the measurement of uncertainty in JAK2 IS NGS data

# Description
Repository contains two files for processing the NGS data:
- mutation_context.py
- plot_data.R

# mutation_context.py
This python script takes BAM files and extracts read information in the context of the mutation site (position 1226 from the reference). The script extracts the following information for each read and outputs the results in a single TSV file:
-	Project ID
-	Sample percentage
-	Run number
-	Read ID
-	Read Orientation (forward or reverse)
-	The distance from the start of the read to the mutation site
-	The distance from the hard-clipped start of the read to the mutation site
-	The distance from the soft-clipped start of the read to the mutation site
-	The distance from the mutation site to the soft-clipped end of the read
-	The distance from the mutation site to the hard-clipped end of the read
-	The distance from the mutation site to the end of the read
-	Allele found at the mutation site
-	The number of insertions found within the read
-	The distance from the mutation site to the nearest insertion
-	The number of deletions found within the read
-	The distance from the mutation site to the nearest deletion
-	The number of SNPs found within the read (excluding the mutation site of interest)
-	The total number of variants found within the read (excluding the mutation site of interest)

The script expects the BAM files to adhere to the following pattern:
- Three digit project number DDD
- Sample percentage between first hyphen and "per"
- "rep" followed by repeat number

Example: 093-0per-JAK2-V617F-PUR-PCR-PROD-2017-03-31-PIA-rep4.30.bwa.JAK2.bam

The script expects the following folder structure:

```
NGS-JAK2
│   README.md
│   mutation_context.py
|   plot_data.R
│
└─── data
    └─── input
    |   |   *BAM files for processing*
    │   |   ...
    │
    └─── output
        │   *Processed BAM files*
        │   ...
        └─── tsv_files
            | *Results stored in TSV format*
            | ...
 ```
 
To run, the user should place the BAM files for processing in the data/input folder, navigate to the root directory and run the script using:

```
python mutation_context.py
```

# plot_data.R
This script performs the analysis on the TSV file produced by the previous script. The plot_data.R script calculates the T allele percentage in the samples whilst changing the constraints each time. The script currently adjusts the minimum allowed distance of the mutation site from the start or end of the read. This can be adjusted by the user to explore the effect of adjusting different parameters.

The script then outputs plots for the T allele percentages, the counts of the minor allele (G or T) and the total read count for each sample.
