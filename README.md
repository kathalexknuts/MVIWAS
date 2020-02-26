# MVTWAS Analysis Pipeline

This file outlines the analysis pipeline for the UKBB application in Knutson and Pan (2020). The key steps in our analysis include:

1. Downloading reference data, such as 1000G
2. Downloading, Reformatting, & Clumping+Thresholding UKBB IDP GWAS summary statistics
3. Acquiring IGAP Summary Statistics & extracting overlapping SNPs from UKBB GWAS  
4. Obtain estimated LD matrices for Stage 1 variants
5. Perform TWAS/MV-TWAS 

## Data Acquisition

**UKBB IDP GWAS Summary Statistics**

GWAS summary statistics on 3,144 Imaging Derived Phenotypes (IDPs) using 9,707 participants have been publicly reported by Elliot et al. [1] (paper: https://www.nature.com/articles/s41586-018-0571-7, resource: http://big.stats.ox.ac.uk/about). We first perform univariate TWAS tests using the summary statistics of 1,578 of these IDPs which were deemed heritable using LDScore Regression (Supp. Table 2 in [1]). Linux commands for download of these GWAS can be found at https://www.dropbox.com/s/qhiftre33pi70xs/BIG_summary_stats_files.xls?dl=0. 

For the sake of this example, we describe our analysis pipeline using a single IDP, namely #0119: T1_FIRST_left_hippocampus_volume. To download, use the wget command from the dropbox link above:

```
wget 'https://www.dropbox.com/s/m7yhcx6h6pqmlgl/0019.txt.gz?dl=0' -O 0019.txt.gz
gunzip 0019.txt.gz
```

This file contains columns "MAF", "BETA", "SEBETA", and "PVAL", and must be merged with the SNP/SNP poisition information, downloadable at https://www.dropbox.com/s/6xcofhwbnyre0s5/positions.txt.gz?dl=0&file_subpath=%2Fpositions.txt. We merge these two files and alter some column names to prepare for LD clumping in plink (specifically RSID -> SNP, PVAL -> P, CHROM -> CHR). 

```
paste ./new.positions.txt ./0019.txt > ./IDP0019.txt
sed -e '1s/RSID/SNP/' -e '1s/PVAL/P/' -e '1s/CHROM/CHR/' ./IDP0019.txt > ./IDP0019tmp.txt
```

**1000 Genomes Reference Panel**

For LD clumping and estimation, we use the phase3 1000G reference panel of 503 subjects of European ancestry [3]. These data can be found at https://www.internationalgenome.org/data#download. Straightforward commands for downloading the EUR data in plink file format via R are given below (taken from http://psoerensen.github.io/qgg/articles/1000genome_tutorial.html)

```
download.file(url="https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz",dest="./1000G_Phase3_plinkfiles.tgz")
system("tar -xvzf 1000G_Phase3_plinkfiles.tgz")
```

## Data Processing: LD clumping and thresholding

We use a clumping radius of 1 Mb and R2 cutoff of 0.1. We extract clumped SNPs from the IDP GWAS. Example plink commands are given below for chromosome 21. We perform LD clumping seperately for each chromosome on every IDP. 

```
./plink --bfile ./1000G_EUR_Phase3_plink/1000G.EUR.QC.21 --clump ./IDP0019tmp.txt --chr 21 --clump-p1 1 --clump-p2 1 --clump-r2 0.10 --clump-kb 1000 --out ./IDP0019 
```

We then extract SNPs with p < 5x10^{-5}. 

```
R code here
```

## Data Processing: Merging with IGAP

We obtain AD GWAS summary statistics from Phase 1 of The International Genomics of Alzheimer’s Project (IGAP) [4]. These data can be downloaded at http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php. We extract the SNPs from the clumped UKBB GWAS. In some cases, some of the IDP variants were missing from IGAP data. We exclude these from analysis. Going forward, we refer to the set of remaining SNPs for IDP 0119 as the Stage 1 SNPs-set. 

```
code for merging SNPs
```

As described in [2], the GWAS effect estimates for the remaining variants for each IDP will be used as weights TWAS/MV-TWAS model. 

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
