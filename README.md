# MV-IWAS Analysis Pipeline

This file outlines the analysis pipeline for the UKBB application in "Implicating Causal Brain Imaging Endophenotypes in Alzheimer’s Disease using Multivariable IWAS and GWAS Summary Data" (Knutson, Deng and Pan, 2020). The key steps in our analysis include:

1. Downloading 1000G reference data
2. Downloading, Reformatting, & Clumping+Thresholding UKBB IDP GWAS summary statistics
3. Acquiring IGAP Summary Statistics & extracting overlapping SNPs from UKBB GWAS  
4. Obtain estimated LD matrices for Stage 1 variants
5. Perform IWAS/MV-IWAS 

## Install Package

```
devtools::install_github("hadley/devtools")
require(devtools)
install_github("kathalexknuts/MVIWAS")
require(MVIWAS)
```

## Data Acquisition

**UKBB IDP GWAS Summary Statistics**

GWAS summary statistics on 3,144 Imaging Derived Phenotypes (IDPs) using 9,707 participants have been publicly reported by Elliot et al. [1] (paper: https://www.nature.com/articles/s41586-018-0571-7, resource: http://big.stats.ox.ac.uk/about). We first perform univariate IWAS tests using the summary statistics of 1,578 of these IDPs which were deemed heritable using LDScore Regression (Supp. Table 2 in [1]). Linux commands for download of these GWAS can be found at https://www.dropbox.com/s/qhiftre33pi70xs/BIG_summary_stats_files.xls?dl=0. 

For the sake of this example, we describe our analysis pipeline using a single IDP, namely #0019: T1_FIRST_left_hippocampus_volume. To download, use the wget command from the dropbox link above:

```
system("wget 'https://www.dropbox.com/s/m7yhcx6h6pqmlgl/0019.txt.gz?dl=0' -O 0019.txt.gz")
system("gunzip 0019.txt.gz")
```

This file contains columns "MAF", "BETA", "SEBETA", and "PVAL", and must be merged with the SNP/SNP poisition information, downloadable at https://www.dropbox.com/s/6xcofhwbnyre0s5/positions.txt.gz?dl=0&file_subpath=%2Fpositions.txt. This file will be called new.positions.txt. We merge these two files and alter some column names to prepare for LD clumping in plink (specifically RSID -> SNP, PVAL -> P, CHROM -> CHR). The PVAL (or renamed P) column gives the -log10 p-value, so transform these back as 10^(-P). 

```
system("cat ./0019.txt | awk -F'\t' 'NR==1{print;next}; {$4=10**(-1*$4); print}' > ./0019tmp.txt")
system("paste ./new.positions.txt ./0019tmp.txt > ./IDP0019.txt")
system("sed -e '1s/RSID/SNP/' -e '1s/PVAL/P/' -e '1s/CHROM/CHR/' ./IDP0019.txt > ./IDP0019tmp.txt")
```

**1000 Genomes Reference Panel**

For LD clumping and estimation, we use the phase3 1000G reference panel [3]. Phased genotypes can be found in vcf file format via the hyperlink at https://www.internationalgenome.org/category/vcf/. We convert these data into plink file format (bed, bim, fam) and extract the 503 subjects of european ancestry, IDs which can be found at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel. In the following document, we use the naming convention "1000G.EUR.#" for these files, where # is replaced with the given chromosome.  

The following commands (written in R) can be used for direct downloading and processing these data (per chromosome) below. For simplicity, we only give example commands for a single chromosome (chr 21), but this should be repeated for all chromosomes for full data.

```
system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")

system("plink --vcf ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out chr21")

system("wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel")

populations <- read.table("./integrated_call_samples_v3.20130502.ALL.panel", header = T)
EUR <- as.character(populations[populations$super_pop == "EUR","sample"])
write.table(cbind(EUR, EUR), "./EUR.txt", col.names = F, row.names = F, quote = F)

system("plink --bfile ./chr21 --keep ./EUR.txt --make-bed --out 1000G.EUR.21")

#remove duplicated variants
system("plink --bfile 1000G.EUR.21 --write-snplist --out ./all_snps")

system("cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist")

system("plink --bfile 1000G.EUR.21 --exclude duplicated_snps.snplist --make-bed --out 1000G.EUR.21.DuplicatesRemoved")
```

## Data Processing: LD clumping and thresholding

We use a clumping radius of 1 Mb and R2 cutoff of 0.1. We extract clumped SNPs from the IDP GWAS. Example plink commands are given below for chromosome 21. We perform LD clumping in parallel for chromosome on every IDP. 

```
system("awk -F'\t' 'NR==1{print;next}$1 == 21' ./IDP0019tmp.txt > ./IDP0019_chr21.txt")

system("plink --bfile ./1000G.EUR.21.DuplicatesRemoved --clump ./IDP0019_chr21.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.10 --clump-kb 1000 --out ./IDP0019_chr21")
```

We then extract SNPs with p < 5x10^{-5}. 

```
clumped <- read.table("./IDP0019_chr21.clumped", header = TRUE)
clumped <- clumped[,-ncol(clumped)]
clumped <- na.omit(clumped[clumped$P <= 5*10^(-5),])
write.table(clumped$SNP, "./clumped_rs.txt", row.names = F, col.names = F, quote = F)

system("awk 'FNR==NR {a[$1]; next}; $2 in a' ./clumped_rs.txt ./IDP0019tmp.txt > ./IDP0019_CPT.txt")
IDP_cpt <- read.table("./IDP0019_CPT.txt", header = FALSE, col.names = c("CHR", "SNP", "POS", "REF", "ALT", "MAF", "BETA.X", "SEBETA.X", "P.X"))
```

## Data Processing: Merging with IGAP

We obtain AD GWAS summary statistics from Phase 1 of The International Genomics of Alzheimer’s Project (IGAP) [4]. These data can be downloaded at http://web.pasteur-lille.fr/en/recherche/u744/igap/igap_download.php. Phase 1 includes 17,008 Alzheimer's disease cases and 37,154 controls (total n = 54162). We first linearize all odds ratios using the approach described in [1]. We then extract the SNPs from the clumped UKBB GWAS. In some cases, some of the IDP variants were missing from IGAP data. We exclude these from analysis. Going forward, we refer to the set of remaining SNPs for IDP 0019 as the Stage 1 SNPs-set. 

```
system("awk 'FNR==NR {a[$1]; next}; $3 in a' ./clumped_rs.txt ./IGAP_summary_statistics/IGAP_stage_1.txt > IGAP_overlap.txt")
IGAP <- read.table("./IGAP_overlap.txt", header = F, col.names = c("CHR", "POS", "SNP", "REF", "ALT", "BETA.Y", "SE.Y", "P.Y"))
IDP_cpt <- IDP_cpt[IDP_cpt$SNP %in% IGAP$SNP,]
IDP_IGAP <- merge(IGAP, IDP_cpt, by = c("SNP", "REF", "ALT", "CHR", "POS"))
```

The resulting data frame IDP_IGAP contains summary statistics for IDP 0019 and AD (from IGAP) after clumping and thresholding of the UKBB IDP data. For this specific example (i.e. IDP 0019, chr 21), we have 2 SNPs remaining - rs4819284 and rs4819210. 

We perform this process in parallel for each chromosome for IDP 0019 to yield a full set of genome-wide variants (specifically # for IDP 0019). 

## Estimating LD matrices from a reference panel

We next estimate LD correlations for the Stage 1 SNP-set using plink using the 1000G reference panel. For fewer than 600 variants, the LD matrix can be implemented using the ld_matrix function from the TwoSampleMR package in R. However, we note that ld estimation using TwoSample MR is based on only 502 EUR subjects so, while these 2 approaches yield *very* similar results, the former may be the more appropriate approach.

**plink implementation ( still written as if in R using system() )**
```
write.table(as.character(IDP_IGAP$SNP), "./IDP_IGAP_SNP.txt", col.names = F, row.names = F, quote = F)
system("plink --bfile ./1000G.EUR.21.DuplicatesRemoved --extract ./IDP_IGAP_SNP.txt --r2 square --write-snplist --out ./IDP0019_chr21")
ld21 <- as.matrix(read.table("./IDP0019_chr21.ld", header = F)); 
snpslist <- as.character((read.table("./IDP0019_chr21.snplist", header = FALSE))[,1])
dimnames(ld21) <- list(snpslist, snpslist)
```

**R package implementation**
```
devtools::install_github("MRCIEU/TwoSampleMR")
library("TwoSampleMR")
ld21 <- na.omit(abs(ld_matrix(as.character(IDP_IGAP$SNP), with_alleles = FALSE))^2)
```

As described in [2], we use a block diagonal LD matrix (22 blocks by chromosome) in IWAS/MV-IWAS. Given a list of the 22 LD matrices for IDP 0119, we use the bdiag function in R from the Matrix package. 

```
library("Matrix")
ld_list <- list(ld1, ld2, ..., ld22)
ZTZ <- as.matrix(bdiag(ld_list))
dimnames(ZTZ) <- list(unlist(lapply(ld_list, rownames)), unlist(lapply(ld_list, rownames)))

IDP_IGAP <- IDP_IGAP[IDP_IGAP$SNP %in% rownames(ZTZ),]
```

## Univariate IWAS with summary statistics

We have now obtained all neccessary data for univariate IWAS of IDP 0119. Before performing the following steps, be sure that the SNP order for the ld matrix is the same as the SNP order for the IDP_IGAP df. As previously noted, the total sample size for the AD GWAS is n= 54162 and the LD estimates are based on nR = 503 EUR subjects from 1000G. 

```
n <- 54162
n_case <- 17008
n_control <- 37154
betaZY <- matrix(IDP_IGAP$BETA.Y, ncol = 1)
se_betaZY <- matrix(IDP_IGAP$SE.Y, ncol = 1)
betaZX <- matrix(IDP_IGAP$BETA.X, ncol = 1)
se_betaZX <- matrix(IDP_IGAP$SEBETA.X, ncol = 1)
corr_mat <- ld21
trait_type <- "Binary"

res_0019 <- mv_iwas_summ(
  betaZY,
  se_betaZY,
  betaZX,
  se_betaZX,
  corr_mat,
  n,
  trait_type,
  n_case,
  n_control
)

save(res_0019, file = "./0019_chr21only.RData")
```

## Multivariate IWAS with summary statistics

Following univariate testing of all heritable UKBB IDPs, we obtain a set of candidate IDPs with univariate p-value < 0.05. We then perform multivariate IWAS using these candidate IDPs seperately for each brain imaging modality group (functional, diffusion, structural). These modality groups can be inferred using the IDP names listed at https://www.dropbox.com/s/qhiftre33pi70xs/BIG_summary_stats_files.xls?dl=0, along with the table giving counts for each modality type in the supplementary material of Elliot et al. Implementation of MV-IWAS almost directly parallels the univariate example given above. Here, however, betaZX is a p x k matrix of weights, where each vector of IDP weights is represented by a column. Ensure that for each column of betaZX, only variants in the associated IDP's SNP-set have non-zero values from the corresponding GWAS effect estimates. All other cells should be set to zero. 

Results for this model are given in [2] and directly compared to the univariate IWAS results.

## References

1. Elliott, L., et al. (2018). Genome-wide association studies of brain imaging phenotypes in UK Biobank. Nature, 562(7726), 210-216.

2. Knutson, Katherine A., Deng, Yangqing, & Pan, Wei. (2020). Implicating causal brain imaging endophenotypes in Alzheimer’s disease using multivariable IWAS and GWAS summary data. NeuroImage, 223, 117347.

3. “A map of human genome variation from population-scale sequencing,” Nature, vol. 467, no. 7319, p. 1061, 2010.

4. J.-C. Lambert et al., “Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for alzheimer’s disease,” Nature Genetics, vol. 45, pp. 1452––1458, 2013.
