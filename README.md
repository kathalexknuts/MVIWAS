# MVTWAS Analysis Pipeline

This file outlines the analysis pipeline for the UKBB application in Knutson and Pan (2020). The key steps in our analysis include:

1. Downloading UKBB IDP GWAS summary statistics
2. Clumping and Thresholding of UKBB GWAS 
3. Acquiring and pre-processing IGAP Summary Statistics
4. Obtain estimated LD matrices for Stage 1 variants
5. 

## Data Acquisition

**UKBB IDP GWAS Summary Statistics**

GWAS summary statistics on 3,144 Imaging Derived Phenotypes (IDPs) using 9,707 participants have been publicly reported by Elliot et al. [1] (paper: https://www.nature.com/articles/s41586-018-0571-7, resource: http://big.stats.ox.ac.uk/about). We first perform univariate TWAS tests using the summary statistics of 1,578 of these IDPs which were deemed heritable using LDScore Regression (Supp. Table 2 in [1]). Linux commands for download of these GWAS can be found at https://www.dropbox.com/s/qhiftre33pi70xs/BIG_summary_stats_files.xls?dl=0. 

For the sake of this example, we describe our analysis pipeline using 1 IDP, INSERT HERE. To download, use the command:

```
insert wget command here
```

## Data Processing: LD clumping and thresholding

We perform LD clumping on each of the UKBB IDP summary statistics using the 1000G reference panel of 503 subjects of European ancestry [3]. These data can be found at https://www.internationalgenome.org/data#download (CHECK THIS HERE!).

For ld clumping, we use a clumping radius of 1 Mb and R2 cutoff of 0.1. We extract clumped SNPs from the IDP GWAS. Example plink commands are given below.

```
plink LD clumping command
```

We load the resulting clumped data in R and extract SNPs with p < 5x10^{-5}. 

```
R code here
```

As described in [2], the GWAS effect estimates for the remaining variants for each IDP will be used as weights TWAS/MV-TWAS model. 

## Data Processing: Merging with IGAP

We obtain AD GWAS summary statistics from Phase 1 of The International Genomics of Alzheimer’s Project (IGAP) [4]. These data can be downloaded at http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php. We extract the SNPs from the clumped UKBB GWAS. In some cases, some of the IDP variants were missing from IGAP data. We exclude these from analysis. Going forward, we refer to the set of remaining SNPs for IDP INSERT NUMBER as the Stage 1 SNPs-set. 

```
code for loading in data and extracting SNPs
```


## Estimating LD matrices from a reference panel

We next estimate LD correlations for the Stage 1 SNP-set. This can be done in plink (given 1000G reference panel) or, for fewer than INSERT NUMBER variants, can be implemented using the INSERT FUNCTION AND PACKAGE NAME HERE. We have confirmed that these two approaches give the same resulting LD matrix.

**plink implementation**
```
plink code
```
**R implementation**
```
R code
```

## Univariate TWAS with summary statistics

We have now obtained all neccessary data for univariate TWAS of IDP INSERT NUMBER. 

```
univariate TWAS code here.
```

The univariate TWAS results for each of these IDPs in given in the Supplementary Materials of [2]. 

## Multivariate TWAS with summary statistics

Following univariate testing of all heritable UKBB IDPs, we obtain the set of NUMBER candidate IDPs with univariate p-value < 0.1. We then perform a multivariate TWAS using these candidate IDPs. Example commands are given below. Here, W is a matrix of weights, where each IDP represents a column. The first command given below ensures that, for each IDP column, only variants in the given IDPs SNP-set have non-zero values from the corresponding GWAS effect estimates. All other cells should be set to zero.  

```
Insert R code here.
```

Results for this model are given in [2] and directly compared to the univariate TWAS results. We discovered NUMBER new IDPs with putative causal effects on AD using the MV-TWAS model over it's univariate analogue.  

## References

1. Elliott, L., et al. (2018). Genome-wide association studies of brain imaging phenotypes in UK Biobank. Nature, 562(7726), 210-216.

2. Knutson, K., Pan, W. (TBD) 

3. “A map of human genome variation from population-scale sequencing,” Nature, vol. 467, no. 7319, p. 1061, 2010.

4. J.-C. Lambert et al., “Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for alzheimer’s disease,” Nature Genetics, vol. 45, pp. 1452––1458, 2013.
