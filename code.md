## Mapping and Variant Calling
```bash
nf-core/sarek

module load Nextflow/25.10.2

cd /path/to/directory/Cr_all_sarek

nextflow pull nf-core/sarek -r 3.8.1

nextflow run nf-core/sarek -r 3.8.1 -profile apptainer -work-dir /path/to/directory/Cr_all_sarek -params-file Cr_params.json -c Cr_custom.config --length-required --with-timeline
# option can be used if run needs to resume after timing out, etc. " -resume "
```

## Variant Filtering
```bash
# Assigning SNP identifiers
bcftools annotate --set-id '%CHROM\_%POS' joint_germline.vcf.gz --output joint_germline_uniqID.vcf.gz --output-type z

vcftools --gzvcf joint_germline_uniqID.vcf.gz --not-chr Crod1.0_mt --recode --recode-INFO-all --out joint_germline_uniqID_noMT

# Preliminary Filtering from sarek output
vcftools --vcf joint_germline_uniqID_noMt.recode.vcf --max-missing 0.80 --minQ 50 --min-meanDP 4 --recode --recode-INFO-all --out joint_mis8q50dp4_uniqID

# High individual missingness
vcftools --vcf joint_mis8q50dp4_uniqID.recode.vcf --missing-indv --out  joint_mis8q50dp4_uniqID

vcftools --vcf  joint_mis8q50dp4_uniqID.recode.vcf --remove Cr_rm_indv_list.txt --recode --recode-INFO-all --out joint_mis8q50dp4_uniqID_rmIndv

# Further Filtering
vcftools --vcf joint_mis8q50dp4_uniqID_rmIndv.recode.vcf --max-missing 0.95 --minQ 200 --min-meanDP 4 --recode --recode-INFO-all --out joint_mis95q200dp4_uniqID_rmIndv

bcftools view --output Cr_mis95q200Adp6_uniqID_rmIndv.vcf.gz --output-type z --exclude-types indels --min-alleles 2 --max-alleles 2 --min-af 0.01:minor --include 'SMPL_MIN(FORMAT/AD) >=6 && SMPL_MAX(FORMAT/AD) <=100' joint_mis95q200dp4_uniqID_rmIndv.recode.vcf

# Identifying SNP loci with high missingness
vcftools --gzvcf Cr_mis95q200Adp6_uniqID_rmIndv.vcf.gz --missing-indv --out  Cr_mis95q200Adp6_uniqID_rmIndv
vcftools --gzvcf Cr_mis95q200Adp6_uniqID_rmIndv.vcf.gz --site-depth --out Cr_mis95q200Adp6_uniqID_rmIndv
vcftools --gzvcf Cr_mis95q200Adp6_uniqID_rmIndv.vcf.gz --site-mean-depth --out Cr_mis95q200Adp6_uniqID_rmIndv
vcftools --gzvcf Cr_mis95q200Adp6_uniqID_rmIndv.vcf.gz--missing-site --out Cr_mis95q200Adp6_uniqID_rmIndv

vcftools --gzvcf Cr_mis95q200Adp6_uniqID_rmIndv --exclude rm_SNPs_list.txt --recode --recode-INFO-all --out Cr_mis95q200Adp6_uniqID_rmIndv_exSNP

# Filtering more 3% missingness
#vcftools --vcf Cr_mis95q200Adp6_uniqID_rmIndv_exSNP.recode.vcf --max-missing 0.97 --recode --recode-INFO-all --out Cr_q200Adp6mis97_uniqID_rmIndv_exSNP

# Identified related individuals and outlier
vcftools --vcf Cr_q200Adp6mis97_uniqID_rmIndv_exSNP.recode.vcf --relatedness --out Cr_q200Adp6mis97_uniqID_rmIndv_exSNP
vcftools --gzvcf Cr_q200Adp6mis97_uniqID_rmIndv_exSNP.recode.vcf --remove rm_related_S194.txt --recode --recode-INFO-all --out Cr_filt_uniqID_exSNP_noRel

# Final File
Cr_filt_uniqID_exSNP_noRel.vcf.gz

# Split by meta-population for subsequent analyses
vcftools --gzvcf Cr_filt_uniqID_exSNP_noRel.vcf.gz --remove ../pop_files/Cr_RA.txt --recode --recode-INFO-all --out Cr_OZNZ
vcftools --gzvcf Cr_filt_uniqID_exSNP_noRel.vcf.gz --remove ../pop_files/Cr_OZNZ.txt --recode --recode-INFO-all --out Cr_RA
vcftools --vcf Cr_OZNZ.recode.vcf --remove ../pop_files/Cr_NZ.txt --recode --recode-INFO-all --out Cr_OZ
vcftools --vcf Cr_OZNZ.recode.vcf --remove ../pop_files/Cr_OZ.txt --recode --recode-INFO-all --out Cr_NZ

```
## LD Pruning and Thinning
```bash

# PopLDdecay 
cd /path/to/Cr_pop_analysis/LD_decay/
#dir="/path/to/Cr_pop_analysis/LD_decay/"

#./PopLDdecay/bin/PopLDdecay -InVCF $dir/../noLDnoHWD/Cr_filt_uniqID_exSNP_noRel.vcf.gz -MaxDist 250 -OutStat $dir/Cr_all.stat.gz
#./PopLDdecay/bin/PopLDdecay -InVCF $dir/Cr_OZ.recode.vcf -MaxDist 250 -OutStat $dir/Cr_OZ.stat.gz
#./PopLDdecay/bin/PopLDdecay -InVCF $dir/Cr_NZ.recode.vcf -MaxDist 250 -OutStat $dir/Cr_NZ.stat.gz
#./PopLDdecay/bin/PopLDdecay -InVCF $dir/Cr_RA.recode.vcf -MaxDist 250 -OutStat $dir/Cr_RA.stat.gz

perl ./PopLDdecay/bin/Plot_MultiPop.pl -inList Cr_multipop.list -output Cr_LD_pops_out -keepR


# LD pruning
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 PLINK/2.00a6.9 BCFtools/1.22-GCC-12.3.0 Python/3.11.6-foss-2023a

cd /path/to/Cr_pop_analysis/LD
plink2 --vcf ../noLDnoHWD/Cr_filt_uniqID_exSNP_noRel.vcf.gz --indep-pairwise 500kb 0.2 --out Cr_LD_all_sites
bcftools view -i 'ID=@Cr_LD_all_sites.prune.in' ../noLDnoHWD/Cr_filt_uniqID_exSNP_noRel.vcf.gz --output Cr_LD_all.vcf.gz --output-type z

cd /path/to/Cr_pop_analysis/LDprune
# Split by populations LD prune dataset
vcftools --gzvcf Cr_LD_all.vcf.gz --remove ../pop_files/Cr_RA.txt --recode --recode-INFO-all --out Cr_LD_OZNZ
vcftools --gzvcf Cr_LD_all.vcf.gz --remove ../pop_files/Cr_OZNZ.txt --recode --recode-INFO-all --out Cr_LD_RA
vcftools --vcf Cr_LD_OZNZ.recode.vcf --remove ../pop_files/Cr_NZ.txt --recode --recode-INFO-all --out Cr_LD_OZ
vcftools --vcf Cr_LD_OZNZ.recode.vcf --remove ../pop_files/Cr_OZ.txt --recode --recode-INFO-all --out Cr_LD_NZ

# Thin dataset 
vcftools --gzvcf Cr_filt_uniqID_exSNP_noRel.vcf.gz --thin 250 --recode --recode-INFO-all --out Cr_thin_all

# Split by populations Thin dataset
vcftools --vcf Cr_thin_all.recode.vcf --remove ../pop_files/Cr_RA.txt --recode --recode-INFO-all --out Cr_thin_OZNZ
vcftools --vcf Cr_thin_all.recode.vcf --remove ../pop_files/Cr_OZNZ.txt --recode --recode-INFO-all --out Cr_thin_RA
vcftools --vcf Cr_thin_OZNZ.recode.vcf --remove ../pop_files/Cr_NZ.txt --recode --recode-INFO-all --out Cr_thin_OZ
vcftools --vcf Cr_thin_OZNZ.recode.vcf --remove ../pop_files/Cr_OZ.txt --recode --recode-INFO-all --out  Cr_thin_NZ
```

## Population Genomic Diversity
```bash

# bash slurm script to run R script
module load R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1
cd /path/to/Cr_pop_analysis/scripts/R_script
srun Rscript ./Cr_hierfstat_slurm.R

# R script
library(vcfR)
library(hierfstat)
library(adegenet)
library(gaston)
setwd("/path/to/Cr_pop_analysis/LD")

Cr_all_data <- "Cr_LD_all.vcf"
Cr_all_vcf <- read.vcfR(Cr_all_data)
Cr_all_gt <- extract.gt(Cr_all_vcf)
Cr_all_gt[Cr_all_gt %in% c("0|0", "0/0")] <- 11
Cr_all_gt[Cr_all_gt %in% c("1/0", "0/1", "1|0", "0|1")] <- 12
Cr_all_gt[Cr_all_gt %in% c("1|1", "1/1")] <- 22
Cr_all_gt <- apply(Cr_all_gt, 2, as.numeric)
pops<-read.table("pops.txt",as.is=T)
Cr_all_gn <- as.data.frame(cbind(pop = as.numeric(as.factor(pops$V1[1:285])), t(Cr_all_gt)))
colnames(Cr_all_gn)[-1] <- paste0("L", 1:nrow(Cr_all_gt))
#saveRDS(Cr_all_gn, "Cr_all_gn.rds")
#write_tsv(Cr_all_gn, "Cr_all_gn.tsv")
#Ht_all_gn <- readRDS("Ht_all_gn.rds")

Cr_all_stats <- basic.stats(Cr_all_gn)
Cr_all_stats$overall
mean_Ho_pop <- colMeans(Cr_all_stats$Ho, na.rm = TRUE)
mean_Hs_pop <- colMeans(Cr_all_stats$Hs, na.rm = TRUE)
mean_Fis_pop <- colMeans(Cr_all_stats$Fis, na.rm = TRUE)
Cr_betas <- betas(Cr_all_gn)
Cr_all_richness <- allelic.richness(Cr_all_gn)
write.csv(Cr_all_richness, "Cr_all_richness.csv")
save.image("Cr_all_hierfstat.RData")

## Visualisation of hierfstat results

setwd("/path/to/Cr_pop_analysis/LD")
library(ggplot2)
library(tidyverse)
library(patchwork)

# Read in table created from hierfstat outputs
data<- read_table("Cr_hierfstat_results.txt")
df <- data.frame(
  pop =data$pop,
  localfst =data$Fst,
  het =data$Ho,
  fis =data$Fis,
  alrich =data$allelicRichness)

fst_plot<- ggplot(df, aes(x = factor(pop) , y = "Local FST", fill = localfst)) +
  geom_tile() +
  scale_fill_gradient2(
    name = NULL,
    low = "plum1",
    mid = "mediumorchid1",
    high = "mediumpurple4",
    midpoint = 0.0) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(family = "serif"),
                             axis.ticks.x = element_blank(), panel.background = element_blank(), plot.background = element_blank())
plot(fst_plot)

het_plot<- ggplot(df, aes(x = factor(pop) , y = "Heterozygosity", fill = het)) +
  geom_tile() +
  scale_fill_gradient2(
    name = NULL,
    low = "plum1",
    mid = "mediumorchid1",
    high = "mediumpurple4",
    midpoint = 0.055) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(family = "serif"),
        axis.ticks.x = element_blank(), panel.background = element_blank(), plot.background = element_blank())
plot(het_plot)

fis_plot<- ggplot(df, aes(x = factor(pop) , y = "FIS", fill = fis)) +
  geom_tile() +
  scale_fill_gradient2(
    name = NULL,
    low = "plum1",
    mid = "mediumorchid1",
    high = "mediumpurple4",
    midpoint = 0.15) +
  labs(x = "Population", y = "") +
  theme(axis.text.x = element_text(size = 14), text = element_text(family = "serif"))  +
  theme(panel.background = element_blank(), plot.background = element_blank())
plot(fis_plot)

alrich_plot<- ggplot(df, aes(x = factor(pop) , y = "Allelic Richness", fill = alrich)) +
  geom_tile() +
  scale_fill_gradient2(
    name = NULL,
    low = "plum1",
    mid = "mediumorchid1",
    high = "mediumpurple4",
    midpoint = 1.26) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(), text = element_text(family = "serif"),
        axis.ticks.x = element_blank(), panel.background = element_blank(), plot.background = element_blank())
plot(alrich_plot)

(alrich_plot/het_plot/fst_plot/fis_plot) +
  plot_layout(heights = c(1, 1)) +  plot_layout(ncol = 1, nrow = 4) +
  theme(plot.margin = margin())


# Private Alleles
# bash script
module load R/4.3.1-gimkl-2022a R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1
cd /path/to/Cr_pop_analysis/scripts/R_script
srun Rscript ./Cr_private_alleles.R

# R script
setwd("/path/to/Cr_pop_analysis/LD")
vcf3_vcfR <- read.vcfR("Cr_LD_all.vcf")
vcf3_gi <- vcfR2genind(vcf3_vcfR)
pop_info <- read.table("pops.txt")
adegenet::pop(vcf3_gi) <- pop_info$V1

## private alleles
PA <- private_alleles(vcf3_gi)
rowSums(PA)
write.csv(PA, "Cr_PA.csv")

# Tajima's D

# Split each sample location examples 
vcftools --vcf Cr_OZ.recode.vcf --keep pop1.txt --recode --recode-INFO-all --out Cr_OZ_1
vcftools --vcf Cr_NZ.recode.vcf --keep pop15.txt --recode --recode-INFO-all --out Cr_NZ_15
vcftools --vcf Cr_RA.recode.vcf --keep pop21.txt --recode --recode-INFO-all --out Cr_RA_21

# Each sample location run in a loop for Tajima's D analysis
module load VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1 PLINK/2.00a6.9 BCFtools/1.22-GCC-12.3.0 Python/3.11.6-foss-2023a
cd /path/to/Cr_pop_analysis/noLD

for file in *.recode.vcf; do
     #Define output filename based on input filename
   output="${file%.recode.vcf}"
 vcftools --vcf "$file" --TajimaD 10000000 --out "$output"_10M
done

## Visualisation of Tajima's D
setwd("/path/to/Cr_pop_analysis/noLD/tajD")
library(ggplot2)
library(tidyverse)
data<- read_table("Cr_allTajD.txt")
df <- data.frame(
  Pops =data$POP,
  Chroms =data$CHROM,
  TajD = data$D)

ggplot(df, aes(x = Pops, y = Chroms, fill = TajD)) +
  geom_tile(color = "black", linewidth = 0.2) +  # thin black borders
  scale_fill_gradient2(  low = "plum1",
                         mid = "mediumorchid1",
                         high = "mediumpurple4", midpoint = -1.95) +
  scale_x_discrete(limits = unique(df$Pops)) +  
  scale_y_discrete(limits = unique(df$Chroms)) +
  labs(x = "Population", y = "Chromosomes") +
  theme_minimal() +
  theme(panel.grid = element_blank())
```

## Population Genomic Structure
```bash

# PCA
module load PLINK/2.00a6.9
cd /path/to/LD
# make bed files
plink2 --vcf Cr_LD_all.vcf --make-bed --out Cr_LD_all
plink2 --vcf Cr_LD_OZ.recode.vcf --make-bed --out Cr_LD_OZ
plink2 --vcf Cr_LD_NZ.recode.vcf --make-bed --out Cr_LD_NZ
plink2 --vcf Cr_LD_RA.recode.vcf --make-bed --out Cr_LD_RA

# run PCA
plink2 --bfile Cr_LD_all --pca --out ./pca/Cr_LD_all
plink2 --bfile Cr_LD_OZ --pca --out ./pca/Cr_LD_OZ
plink2 --bfile Cr_LD_NZ --pca --out ./pca/Cr_LD_NZ
plink2 --bfile Cr_LD_RA --freq --out Cr_LD_RA
plink2 --bfile Cr_LD_RA --read-freq Cr_LD_RA.afreq --pca --out ./pca/Cr_LD_RA

## Visualisation of PCA RStudio
library(tidyverse)
library(RcppCNPy)
library(ggplot2)
library(ggrepel)
setwd("/path/to/Cr_pop_analysis/LD/pca")

### All
pca <- read_table2("Cr_LD_all.eigenvec", col_names = TRUE)
eigenval <- scan("Cr_LD_all.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
pca <- pca %>% relocate(IID, .after = last_col())
pops <- read.table("Cr_all_pops.txt", header = TRUE)
pca <- as_tibble(data.frame(pca, pops))
pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
breaks <- c("ByronBay",	"SWFishRocks",	"Sydney",	"SydneyHistoric",	"JervisBay",	"CapeHowe",
"DealIsland",	"DealIslandHistoric",	"Gardens",	"ThumbsHistoric",	"Thumbs",	"LordHoweIsland",	
"LordHoweIslandHistoric",	"NorfolkIsland",	"NorthCape",	"CavalliIslands",	"PoorKnightsIslands",	
"MokohinauIslands",	"MercuryIslands",	"Whakaari",	"L'EsperanceRock",	"RangitāhuIsland")
highlight <- c("ByronBay","SWFishRocks", "SydneyHistoric","JervisBay", "CapeHowe", "DealIslandHistoric", "ThumbsHistoric", "LordHoweIslandHistoric")
fill_cols <- c(
  "ByronBay"="#721E17", "SWFishRocks"="#95211B", "Sydney"="#B8221E", "SydneyHistoric"="#DA2222", "JervisBay"="#DF4828", "CapeHowe"="#E4632D",
  "DealIsland"="#E67932", "DealIslandHistoric"="#E78C35", "Gardens"="orange", "ThumbsHistoric"="goldenrod", "Thumbs"="goldenrod1", "LordHoweIsland"="gold2",
  "LordHoweIslandHistoric"="gold", "NorfolkIsland"="yellow2", "NorthCape"="#549EB3", "CavalliIslands"="#59A5A9", "PoorKnightsIslands"="#60AB9E",
  "MokohinauIslands"="#69B190", "MercuryIslands"="skyblue",  "Whakaari"="dodgerblue", "L'EsperanceRock"="purple4",  "RangitāhuIsland"="slateblue3")
outline_cols <- fill_cols
outline_cols[highlight] <- "black"
b <- ggplot(pca, aes(x = PC1, y = PC2,fill = pop, colour = pop)) + geom_point(shape = 21, size = 3, stroke = 0.8) 
b <- b + scale_fill_manual(breaks = breaks, values = fill_cols) + scale_colour_manual(breaks = breaks, values = outline_cols)
b<- b  +  theme_light() +  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 18), 
                                 axis.text.y = element_text(size = 18), legend.title = element_text(size = 0),
                                 legend.text = element_text(size = 18), text = element_text(family = "serif"), axis.title = element_text(size = 18)) 
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

### Aotearoa
pca <- read_table2("Cr_LD_NZ.eigenvec", col_names = TRUE)
eigenval <- scan("Cr_LD_NZ.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
pca <- pca %>% relocate(IID, .after = last_col())
pops <- read.table("Cr_NZ_pops.txt", header = TRUE)
pca <- as_tibble(data.frame(pca, pops))


pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
breaks <- c("NorthCape",	"CavalliIslands",	"PoorKnightsIslands",	"MokohinauIslands",	
              "MercuryIslands",	"Whakaari")
b <- ggplot(pca, aes(PC1, PC2, col = pop))
b <- b + geom_point(size = 3) +  labs(x = "PC1", y = "PC2") + 
  scale_colour_manual(breaks = breaks, 
                      values = c("NorthCape"="#549EB3", "CavalliIslands"="#59A5A9", "PoorKnightsIslands"="#60AB9E", 
                                 "MokohinauIslands"="#69B190", "MercuryIslands"="skyblue", "Whakaari"="dodgerblue")) + 
  theme(panel.grid = element_blank())
b <- b + theme_light() + theme(panel.grid = element_blank())
b<- b  +  theme_light() +  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 18), 
                                 axis.text.y = element_text(size = 18), legend.title = element_text(size = 0),
                                 legend.text = element_text(size = 18), text = element_text(family = "serif"), axis.title = element_text(size = 18)) 
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

### Australia
pca <- read_table2("Cr_LD_OZ.eigenvec", col_names = TRUE)
eigenval <- scan("Cr_LD_OZ.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
pca <- pca %>% relocate(IID, .after = last_col())
pops <- read.table("Cr_OZ_pops.txt", header = TRUE)
pca <- as_tibble(data.frame(pca, pops))


pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
breaks <- c("ByronBay",	"SWFishRocks",	"Sydney",	"SydneyHistoric",	"JervisBay",	"CapeHowe",
            "DealIsland",	"DealIslandH",	"Gardens",	"ThumbsHistoric",	"Thumbs",	"LordHoweIsland",	
            "LordHoweIslandHistoric",	"NorfolkIsland")
highlight <- c("ByronBay","SWFishRocks", "SydneyH","JervisBay", "CapeHowe", "DealIslandH", "ThumbsH", "LordHoweH")
fill_cols <- c(
  "ByronBay"="#721E17", "SWFishRocks"="#95211B", "Sydney"="#B8221E", "SydneyHistoric"="#DA2222", "JervisBay"="#DF4828", "CapeHowe"="#E4632D",
  "DealIsland"="#E67932", "DealIslandHistoric"="#E78C35", "Gardens"="orange", "ThumbsHistoric"="goldenrod", "Thumbs"="goldenrod1", "LordHoweIsland"="gold2",
  "LordHoweIslandHistoric"="gold", "NorfolkIsland"="yellow2")
outline_cols <- fill_cols
outline_cols[highlight] <- "black"
b <- ggplot(pca, aes(x = PC1, y = PC2,fill = pop, colour = pop)) + geom_point(shape = 21, size = 3, stroke = 0.8) 
b <- b + scale_fill_manual(breaks = breaks, values = fill_cols) + scale_colour_manual(breaks = breaks, values = outline_cols)
b<- b  +  theme_light() +  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 18), 
                                 axis.text.y = element_text(size = 18), legend.title = element_text(size = 0),
                                 legend.text = element_text(size = 18), text = element_text(family = "serif"), axis.title = element_text(size = 18)) 
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

### Rangitāhua
pca <- read_table2("Cr_LD_RA.eigenvec", col_names = TRUE)
eigenval <- scan("Cr_LD_RA.eigenval")
# sort out the pca data
# remove nuisance column
pca <- pca[,-1]
pca <- pca %>% relocate(IID, .after = last_col())
pops <- read.table("Cr_RA_pops.txt", header = TRUE)
pca <- as_tibble(data.frame(pca, pops))


pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()
cumsum(pve$pve)
b <- ggplot(pca, aes(PC1, PC2, col = pop))
b <- b + geom_point(size = 3) +  labs(x = "PC1", y = "PC2") + 
  scale_colour_manual(breaks = c( "L'EsperanceRock",	"RangitāhuIsland"), 
                      values = c( "L'EsperanceRock"="purple4", "RangitāhuIsland"="slateblue3")) + theme(panel.grid = element_blank())
b<- b  +  theme_light() +  theme(panel.grid = element_blank(), axis.text.x = element_text(size = 18), 
                                 axis.text.y = element_text(size = 18), legend.title = element_text(size = 0),
                                 legend.text = element_text(size = 18), text = element_text(family = "serif"), axis.title = element_text(size = 18)) 
b + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) +  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

# Admixture
# bash script
module load R/4.3.1-gimkl-2022a R-bundle-Bioconductor/3.17-gimkl-2022a-R-4.3.1
cd /path/to/Cr_pop_analysis/scripts/R_script
srun Rscript ./Cr_LEA_admix_slurm.R

# R script
library(LEA)
library(vcfR) 
setwd("/path/to/Cr_pop_analysis/LD/")

### All
project = load.snmfProject("Cr_LD_all.snmfProject")
# plot cross-entropy criterion of all runs of the project
plot(project, cex = 1.2, col = "lightblue", pch = 19)
best = which.min(cross.entropy(project, K = 7))
q <- Q(project, K = 7, run = best)
pop<-read.table("Cr_all_pops.txt",as.is=T)
main_cluster <- apply(q, 1, which.max)
max_prop     <- apply(q, 1, max)
ord <- order(pop[,1], main_cluster, -max_prop)
cols <- c( "orchid4", "aquamarine4", "plum1","mediumpurple4" , "paleturquoise3", "purple","magenta2")
par(family = "serif")
barplot(t(q)[,ord],col= cols,space=0,border=NA,xlab="Individuals",
        ylab="Admixture proportions",  cex.axis = 1.8, font.axis = 1,  cex.lab = 1.8, font.lab = 1) + theme(text = element_text(family = "serif"))
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T, srt = 0, cex=1.8)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=3.0)

### Rangitāhua
output = vcf2geno("Cr_HWD_LD_noOuts_RA.vcf")
RA_project = snmf("Cr_HWD_LD_noOuts_RA.geno",
K = 1:10, 
entropy = TRUE, 
repetitions = 5,
project = "new",
CPU = 8, seed = 7)
project = load.snmfProject("Cr_LD_RA.recode.snmfProject")
plot(project, cex = 1.2, col = "lightblue", pch = 19)
best = which.min(cross.entropy(project, K = 2))
q <- Q(project, K = 2, run = best)
pop<-read.table("Cr_RA_pops.txt",as.is=T)  
main_cluster <- apply(q, 1, which.max)
max_prop     <- apply(q, 1, max)
ord <- order(pop[,1], main_cluster, -max_prop)
cols <- c("orchid4", "purple", "mediumpurple4", "plum1", "orchid4", "magenta2")
par(family = "serif")
barplot(t(q)[,ord],col= cols,space=0,border=NA,xlab="Individuals",
        ylab="Admixture proportions",  cex.axis = 1.8, font.axis = 1,  cex.lab = 1.8, font.lab = 1) + theme(text = element_text(family = "serif"))
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T, srt = 0, cex=1.8)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=3.0)
dev.off()

### Aotearoa
output = vcf2geno("Cr_HWD_LD_noOuts_NZ.vcf")
NZ_project = snmf("Cr_HWD_LD_noOuts_NZ.geno",
                 K = 2:10, 
                entropy = TRUE, 
                repetitions = 10,
                project = "new",
                CPU = 8, seed = 30)
project = load.snmfProject("Cr_LD_NZ.recode.snmfProject")
plot(project, cex = 1.2, col = "lightblue", pch = 19)
best = which.min(cross.entropy(project, K = 4))
q <- Q(project, K = 4, run = best)
pop<-read.table("Cr_NZ_pops.txt",as.is=T)    
main_cluster <- apply(q, 1, which.max)
max_prop     <- apply(q, 1, max)
ord <- order(pop[,1], main_cluster, -max_prop)
cols <- c( "purple", "mediumpurple4", "plum1", "orchid4")
par(family = "serif")
barplot(t(q)[,ord],col= cols,space=0,border=NA,xlab="Individuals",
        ylab="Admixture proportions",  cex.axis = 1.8, font.axis = 1,  cex.lab = 1.8, font.lab = 1) + theme(text = element_text(family = "serif"))
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T, srt = 0, cex=1.8)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=3.0)
dev.off()

### Australia
output = vcf2geno("Cr_HWD_LD_noOuts_OZ.vcf")
OZ_project = snmf("Cr_HWD_LD_noOuts_OZ.geno",
                 K = 2:22, 
                entropy = TRUE, 
               repetitions = 5,
              project = "new",
             CPU = 8, seed = 30)
project = load.snmfProject("Cr_LD_OZ.recode.snmfProject")
plot(project, cex = 1.2, col = "lightblue", pch = 19)
best = which.min(cross.entropy(project, K = 4))
q <- Q(project, K = 4, run = best)
pop<-read.table("Cr_OZ_pops.txt",as.is=T)
main_cluster <- apply(q, 1, which.max)
max_prop     <- apply(q, 1, max)
ord <- order(pop[,1], main_cluster, -max_prop)
cols <- c( "purple" ,"orchid4" ,"plum1" , "mediumpurple4")
cols <- c( "paleturquoise3", "mediumpurple4", "plum1", "orchid4", "aquamarine4", "purple","seagreen3","magenta3")
par(family = "serif")
barplot(t(q)[,ord],col= cols,space=0,border=NA,xlab="Individuals",
        ylab="Admixture proportions",  cex.axis = 1.8, font.axis = 1,  cex.lab = 1.8, font.lab = 1) + theme(text = element_text(family = "serif"))
text(tapply(1:nrow(pop),pop[ord,1],mean),-0.05,unique(pop[ord,1]),xpd=T, srt = 0, cex=1.8)
abline(v=cumsum(sapply(unique(pop[ord,1]),function(x){sum(pop[ord,1]==x)})),col=1,lwd=3.0)
dev.off()

# Isolation by Distance

# First run Identity by State
# create bed files and update "family ID" by sample location
# bash script
module load PLINK/2.00a6.9 
cd /path/to/Cr_pop_analysis/LD
plink2 --vcf Cr_LD_all.vcf.gz --make-bed --out Cr_LD_all
plink2 --bfile Cr_LD_all --update-ids update_fam.txt --make-bed --out Cr_LD_all_upFam
plink2 --vcf Cr_LD_NZ.recode.vcf --make-bed --out Cr_LD_NZ
plink2 --vcf Cr_LD_OZ.recode.vcf --make-bed --out Cr_LD_OZ
plink2 --vcf Cr_LD_RA.recode.vcf --make-bed --out Cr_LD_RA
plink2 --bfile Cr_LD_NZ --update-ids update_fam.txt --make-bed --out Cr_LD_NZ_upFam
plink2 --bfile Cr_LD_OZ --update-ids update_fam.txt --make-bed --out Cr_LD_OZ_upFam
plink2 --bfile Cr_LD_RA --update-ids update_fam.txt --make-bed --out Cr_LD_RA_upFam

module load PLINK/1.09b6.16
plink --bfile Cr_LD_all_upFam --distance square 1-ibs --allow-extra-chr --out Cr_LD_all_upFam
plink --bfile Cr_LD_NZ_upFam --distance square 1-ibs --allow-extra-chr --out Cr_LD_NZ_upFam
plink --bfile Cr_LD_OZ_upFam --distance square 1-ibs --allow-extra-chr --out Cr_LD_OZ_upFam
plink --bfile Cr_LD_RA_upFam --distance square 1-ibs --allow-extra-chr --out Cr_LD_RA_upFam

# create IBS means
# R scrip
setwd("/path/to/Cr_pop_analysis/LD/IBS")
library(tidyverse)

# Load data
matrix_data <- read.table("Cr_LD_all_upFam.mdist", header = FALSE)
ids <- read.table("Cr_LD_all_upFam.mdist.id", header = FALSE, col.names = c("FID", "IID"))
pop_assigned <- read.table("pop_assign.txt", header = FALSE, col.names = c("FID", "IID", "Population"))

# Assign population labels to rows and columns
names(matrix_data) <- pop_assigned$Population
matrix_data$Pop_Row <- pop_assigned$Population

# Convert to long format and calculate means
ibs_means <- matrix_data %>%
  pivot_longer(cols = -Pop_Row, names_to = "Pop_Col", values_to = "IBS") %>%
  group_by(Pop_Row, Pop_Col) %>%
  summarize(Mean_IBS = mean(IBS), .groups = 'drop')

print(ibs_means)
write.table(ibs_means, file = "Cr_LD_IBS_means.txt")

# Run marmap to get geographic distance
module load R-Geo/4.3.2-foss-2023a
cd /path/to/Cr_pop_analysis/Cr_IBD
srun Rscript ./Cr_IBDistance.R

# R script 
setwd("/nesi/nobackup/ga03714/Cr_raw_lcWGS_data/Cr_pop_analysis/Cr_IBD")
library(ade4)
library(ggplot2)
library(marmap)
library(ncdf4)
Cr_tas.sites <- read.table("Cr_sample_locations", header = TRUE)
tasman<- readGEBCO.bathy("gebco_2026_n-27.0_s-45.0_w145.0_e185.0.nc")
trans<- trans.mat(tasman)
Cr_out.tas.dist <- lc.dist(trans, Cr_tas.sites, res = "dist")
Cr_out.tas.path <- lc.dist(trans, Cr_tas.sites, res = "path")
save.image("/path/to/Cr_pop_analysis/Cr_IBD/Cr_IBD_dist.RData")

load(file = "/path/to/Cr_pop_analysis/Cr_IBD/Cr_IBD_dist.RData")
plot(tasman, n = 1, image = TRUE, main = "Tasman Region")
points(Cr_tas.sites, pch = 21, col = "purple", bg = col2alpha("purple", .9), cex = 1.2)
lapply(Cr_out.tas.path, lines, col = "purple", lwd = 0.4, lty = 2)-> dummy
Cr_dist_matrix = as.matrix(Cr_out.tas.dist)
write.table(Cr_dist_matrix, file ="Cr_dist_matrix.txt")
view(Cr_dist_matrix.txt)


# IBD analysis

# IBD with FST
library(ade4)
library(ggplot2)
# Load matrices
fst_raw <- read.table("Cr_fst_matrix.txt",  row.names = 1, check.names = FALSE)
geo_raw <- read.table("Cr_dist_matrix.txt", row.names = 1, check.names = FALSE)

# Convert data frames into matrix format
fst_mat <- as.matrix(fst_raw)
geo_mat <- as.matrix(geo_raw)

# Apply standard transformations for 2D habitats (Rousset, 1997)
# Genetic distance: Fst / (1 - Fst)
fst_trans <- fst_mat / (1 - fst_mat)

# Geographic distance: log(distance) 
# use offset for exact zero/same pop or change 0 to 1 before log transformation to avoid -Inf values.
geo_mat[geo_mat == 0] <- 1  
geo_log <- log(geo_mat)

#Convert matrices to 'dist' class objects for the Mantel test
gen_dist <- as.dist(fst_trans)
geo_dist <- as.dist(geo_mat)

# Run the Mantel Test
# Using 9999 permutations for robust p-value estimation
set.seed(123) # For reproducibility of permutations
Cr_mantel_FST_result <- mantel.rtest(geo_dist, gen_dist, nrepet = 9999)

# Print the Summary Results to Console
print(Cr_mantel_FST_result)

# Plot
plot(geo_dist, gen_dist, 
     pch = 16, 
     xlab = "Geographic Distance (km)", 
     ylab = "Fst / (1 - Fst)",
     main = "Isolation-by-Distance Plot")
abline(lm(gen_dist ~ geo_dist), col = "purple", lwd = 2)


# IBD with IBS
# Load matrices
ibs_raw <- read.table("Cr_IBS_mean_matrix.txt",  row.names = 1, check.names = FALSE)
geo_raw <- read.table("Cr_dist_matrix.txt", row.names = 1, check.names = FALSE)

# Convert data frames into matrix format
ibs_mat <- as.matrix(ibs_raw)
geo_mat <- as.matrix(geo_raw)

# Geographic distance: log(distance) 
# use offset for exact zero/same pop or change 0 to 1 before log transformation to avoid -Inf values
geo_mat[geo_mat == 0] <- 1  
geo_log <- log(geo_mat)

# Convert matrices to 'dist' class objects for the Mantel test
gen_dist <- as.dist(ibs_mat)
geo_dist <- as.dist(geo_mat)

# Run the Mantel Test
# Using 9999 permutations for robust p-value estimation
set.seed(123) # For reproducibility of permutations
Cr_mantel_IBS_result <- mantel.rtest(geo_dist, gen_dist, nrepet = 9999)

# Print the Summary Results to Console
print(Cr_mantel_IBS_result)

# Plot
plot(geo_dist, gen_dist, 
     pch = 16, 
     xlab = "Geographic Distance", 
     ylab = "Identity by State",
     main = "Isolation-by-Distance Plot")
abline(lm(gen_dist ~ geo_dist), col = "purple", lwd = 2)
```

## Directionality
```bash
#' loads a data set from a snapp and coords file
#'
#'
#' This is a function that reads the old data format. Presented for backwards
#' compatibility
#' @param snapp.file a file in snapp data format
#' @param coords.file a file containing coordinates and other info 
#'      about the sample
#' @param n.snp the maximum number of snp to be read. if -1, all
#'      snps are read
#' @param ploidy the ploidy of the organism, 1 for haploids, 2 for diploid
#' @param ... further arguments passed to load.coord.file
#' @return an object of type origin.data with entries `genotypes`,
#'      containing the genetic data and entry `coord` containing location
#'      data
#' @example examples/example_1.r
#' @export
#' 

# raw.data <- load.data.snapp(snapp.file="GT.snapp",
#                             coords.file="coord.csv",
#                             n.snp=-1,
#                             sep=',', ploidy=2)

library(snpStats)
library(survival)
library(Matrix)
library(geosphere)

load.data.snapp <- function(snapp.file, coords.file, n.snp=-1, 
                            ploidy=2, ...){
  f <- load.snapp.file(snapp.file, n.snp=n.snp)
  coords <- load.coord.file(coords.file, sep=",")
  raw.data <- check.missing(f, coords)
  raw.data <- set.outgroups(raw.data, ploidy)
  class(raw.data) <- 'origin.data'
  return(raw.data)
}

# Add this function for plink input
load.data.plink <- function(plink.file, coords.file, n.snp=NULL,
                            ploidy=2, ...) {
  if (!is.null(n.snp)) {
    f <- read.plink(bed = plink.file, select.snps = 1:n.snp)
  }
  else {
    f <- read.plink(bed = plink.file)
  }
  coords <- load.coord.file(coords.file, sep=",")
  # raw.data <- check.missing(f, coords)
  raw.data <- set.outgroups(raw.data, ploidy)
  # class(raw.data) <- 'origin.data'
  # return(raw.data)
}


#' Reads data from SNAPP file
#' 
#' Deprecated, use plink input format instead
#' @example examples/example_1.r
load.snapp.file <- function(snp.file=snapp.file,
                            n.snp=-1){
  data <- read.table(snp.file, sep =",", header=F, strings=F,
                     na.strings="?", nrow=n.snp, row.names=1)
  ids <- rownames(data)
  # data[1:10, 1:10]
  raw.data <- list()
  data <- apply(as.matrix(data),2,as.numeric)
  rownames(data) <- ids
  #apply(as.matrix(data)[1:10, 1:10],2,as.numeric)
  raw.data$genotypes <- as(data, "SnpMatrix")
  rownames(raw.data$genotypes) <- ids
  raw.data$fam <- rownames(data)
  return(raw.data)
}


#' Loads a coordinate data set
#'
#' Loads a file that specifies location data and other sample specific information
#' A
#'
#' @param file: the name of the file which the data are to be read from.
#' @param ... Additional arguments passed to read.table
#' @return A data.frame with columns sample, longitude, latitude, region, 
#'outgroup and country 
#' @example examples/example_1.r
#' @export

# file <-  "coord.csv"
load.coord.file <- function(file, ...){
  coords <- read.table(file, header=T, strings=F, ...)
  names(coords)[1] <- 'id'
  coords[,1] <- as.character(coords[,1])
  return(coords)
}

#' Cleans genotype and coordinate files for consistency
#'
#' 
#'
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @param snp.data output from load.plink.file
#' @param coord.data output from load.coord.file
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export

# snp.data <- f
# coord.data = coords
check.missing <- function(snp.data, coord.data){
  
  ids <- rownames(snp.data$genotypes)
  to.remove <- !unlist(lapply(ids,
                              function(x)x%in%coord.data$id))
  res.data <- list()
  res.data$coords <- coord.data
  res.data$genotypes <- snp.data$genotypes[!to.remove,]
  ids <- ids[!to.remove]
  
  data.ordering <- res.data$coords$id
  
  res.data$genotypes <- res.data$genotypes[data.ordering,]
  
  return(res.data)
}


#' Sets columns of outgroup individuals
#'
#' As psi requires the derived state of an allele to be known, this function
#' provides an easy way to polarize SNP with outgroup individuals. From the
#' outgroup columns in the coord file, outgroup individuals are removed
#' and SNPs swithced if necessary
#'
#' @param data output from check missing
#' @param ploidy the number of copies per individual
#' @return An object with entries 'genotypes' and 'coords' containing the
#'    input data sets
#'
#' @example examples/example_1.r
#' @export
set.outgroups <- function(data, ploidy=2
){
  if(is.null(data$coords$outgroup)){ 
    return(data)
  }
  
  outgroup.columns <- which(as.logical(data$coords$outgroup))
  n.outgroups <- length(outgroup.columns)
  if( n.outgroups > 0 ){
    data$outgroup <- data$genotypes[outgroup.columns,]
    data$outgroup <- as(data$outgroup, 'numeric')
    data$genotypes <- data$genotypes[-1* outgroup.columns,]
    
    n.alleles <- colSums(!is.na(data$outgroup))
    derived <- colSums(data$outgroup==ploidy, na.rm=T) == n.alleles
    anc <- colSums(data$outgroup==0, na.rm=T)== n.alleles
    no.outgrp <- is.na(derived+anc) | !(anc | derived)  
    
    data$genotypes <- switch.alleles(data$genotypes,derived)
    
    data$genotypes <- data$genotypes[,!no.outgrp ]
    data$coords <- data$coords[-outgroup.columns,]
  }
  
  return(data)
}




#' Sets up the popultion structure data sets
#'
#' Groups up raw data into a data structure that 
#' 
#' @param raw.data Output from set.outgroups or check.data
#' @param ploidy the ploidy of organisms
#' @return A list with entries
#' - data : an n x p matrix with derived allele counts
#' - coords : an p x 3 matrix of coords
#' - n : number of populations
#' @example examples/example_1.r
#' @export
make.pop <- function(raw.data, ploidy=2){
  if(is.null(raw.data$coord$pop)){
    raw.data<- make.pops(raw.data)
  }
  pop <- list()
  pop$data <- make.pop.data(raw.data)[,-1]
  pop$ss <- make.pop.ss(raw.data, ploidy)[,-1]
  pop$coords <- make.pop.coords(raw.data)
  pop$coords$hets <- rowMeans(pop$data/pop$ss, na.rm=T)
  pop$n <- nrow(pop$coords)
  class(pop) <- 'population'
  return(pop)
}
#hets <- (coords$hets - min(coords$hets))/(max(coords$hets) - min(coords$hets))

make.pop.data <- function(raw.data){
  data <- as(raw.data$genotypes, 'numeric')
  aggregate(data, by=list(raw.data$coords$pop), FUN=sum, na.rm=T)
}
make.pop.ss <- function(raw.data, ploidy=2){
  data <- as(raw.data$genotypes, 'numeric')
  aggregate(data, by=list(raw.data$coords$pop), 
            FUN=function(x)ploidy*sum(!is.na(x)))
}

make.pop.coords <- function(raw.data){
  lat <- aggregate(raw.data$coords$latitude, 
                   by=list(raw.data$coords$pop), 
                   FUN=mean)
  names(lat) <- c("pop", "latitude")
  long <- aggregate(raw.data$coords$longitude, 
                    by=list(raw.data$coords$pop), 
                    FUN=mean)
  names(long) <- c("pop", "longitude")
  pop.coords <- merge(lat,long)
  region.id <- match(pop.coords$pop, raw.data$coords$pop)
  pop.coords$region <- raw.data$coords$region[region.id]
  return(pop.coords)
}


#' Puts individuals into populations for further analysis
#'
#' There are two main modes of grouping individuals into populations: 
#' The first one is based on the location, all individuals with the 
#' same sample location are assigned to the same population.
#' the other one is based on a column `population` in 
#'
#'
#' @param raw.data Output from set.outgroups or check.data
#' @param mode One of 'coord' or 'custom'. If 'custom', populations are
#'   generated from the population column in the coords file. If 'coord',
#'   Individuals are grouped according to their coordinates
#' @return A list with each entry being a population
make.pops <- function(raw.data, mode='coord'){
  if(mode == 'coord'){
    s <- c('latitude', 'longitude')
    rowSums <- function(x)200*x[,1]+x[,2]
    raw.data$coords$pop <-  match(rowSums(raw.data$coords[,s]), 
                                  rowSums(raw.data$coords[,s]))
  } 
  return(raw.data)
}

# pop <- pop_sim
# bin= TRUE

flatten <- function(x) {
  if (!inherits(x, "list")) return(list(x))
  else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
}


get.all.psi.mc.bin <- function(pop, n=2,cores=10){
  #this function calculates psi for columns i,j, both resampled down to
  # n samples
  
  n.pops <- pop$n
  
  ## get counts for all pops
  counts <- mclapply(1:n.pops,function(i){
    ni <- unlist(pop$ss[i,])
    fi <- unlist(pop$data[i,])
    list(ni,  fi)
  },mc.cores = cores)
  
  mat = matrix( 0, nrow=n.pops, ncol=n.pops )
  pvals = matrix( 0, nrow=n.pops, ncol=n.pops)
  cnts = matrix( 0, nrow=n.pops, ncol=n.pops)
  
  ## get psi
  # i <- 1
  # j <- 57
  out <- do.call(rbind, flatten(mclapply(1:(n.pops-1),function(i){
    
    lapply((i+1):n.pops, function(j){
      #print(c(i,j))
      c(i,j,get.psi.bin(counts[[i]][[1]], counts[[j]][[1]], counts[[i]][[2]], counts[[j]][[2]]))
      
    })
    
  },mc.cores=cores)))
  
  
  mat[out[,1:2]] <- out[,3]
  mat[lower.tri(mat)] <- -t(mat)[lower.tri(mat)]
  
  
  pvals[out[,1:2]] <- out[,4]
  pvals[lower.tri(pvals)] <- t(pvals)[lower.tri(pvals)]
  diag(pvals) <- NA
  cnts[out[,1:2]] <- out[,5]
  cnts[lower.tri(cnts)] <- t(cnts)[lower.tri(cnts)]
  
  return(list(psi=mat,pvals=pvals,cnts=cnts))
}

#list(psi=NA,pvals=NA,cnts=NA)
get.psi.bin <- function (ni, nj, fi, fj){
  n = 2
  fn <- cbind( fi, ni, fj, nj)
  
  tbl <- table(as.data.frame(fn))
  tbl <- as.data.frame( tbl )
  tbl <- tbl[tbl$Freq > 0,]
  
  tbl$fi <- as.integer(as.character( tbl$fi ))
  tbl$fj <- as.integer(as.character( tbl$fj ))
  tbl$ni <- as.integer(as.character( tbl$ni ))
  tbl$nj <- as.integer(as.character( tbl$nj ))
  
  to.exclude <- tbl$fi == 0 | tbl$fj == 0 | tbl$ni < n | tbl$nj < n
  
  tbl <- tbl[! to.exclude, ]
  
  if(nrow(tbl)==0){ return(NaN)}
  
  poly.mat <- matrix(0,nrow=n+1, ncol=n+1)
  poly.mat[2:(n+1),2:(n+1)] <- 1
  poly.mat[n+1,n+1] <- 0
  
  #psi.mat is the contribution to psi for each entry
  psi.mat <- outer(0:n,0:n,FUN=function(x,y)(y-x))
  psi.mat[1,] <- 0
  psi.mat[,1] <- 0
  
  f.contribution <- function(row, b=n){
    a <- 0:b
    f1 <- row[1]
    n1 <- row[2]
    f2 <- row[3]
    n2 <- row[4]
    cnt <- row[5]
    q1 <- choose(b, a) * choose(n1-b, f1-a)/choose(n1,f1)
    q2 <- choose(b, a) * choose(n2-b, f2-a)/choose(n2,f2)
    return( cnt * outer(q2, q1) )
  }
  
  
  resampled.mat <- matrix(rowSums(apply(tbl, 1,f.contribution)),nrow=n+1)
  sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat)
  
  cnts <- c(as.integer(resampled.mat[2,3]),as.integer(resampled.mat[3,2]))
  
  bin <- ifelse(sum(cnts)>0, binom.test(cnts)$p.value, 1)
  
  return( c(sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) ,bin,sum(cnts)))
}

#' preps data for tdoa analysis
#'
#' makes a data structure of format xi, yi, xj, yj, psi
#' 
#' @param pop.coords coordination data with latitude and longitude
#' @param all.psi matrix with psi values
#' @param region region or set of regions to run analysis on
#' @param countries countries to restrict analysis to
#' @param xlen number of x points
#' @param ylen number of y points
#'
prep.tdoa.data <- function(coords, psi){
  locs <- coords[,c('longitude', 'latitude')]
  n.locs <- nrow(locs)
  
  tdoa.data <- c()
  for(i in 1:n.locs){
    for(j in (i+1):n.locs){
      if( i>=j | j >n.locs) break
      tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
    }
  }
  
  tdoa.data <- matrix(unlist(tdoa.data), ncol=5)
  
  return(tdoa.data)
}

div_from_sfs <- function(pops, GTs,seq_len=1e6){
  rbindlist(lapply(unique(pops),function(pop){
    sfs <- sfs1D(pops,pop,GTs,seq_len)
    n=length(sfs)
    
    # Calculate n-1th harmonic number to use for Whatterson's theta estimation
    numChrom=n-1
    harmonicNumber = 0
    for (j in 1:(numChrom - 1)) {
      harmonicNumber = harmonicNumber + 1.0/j
    }
    
    # creates weights for each class. This assumes the spectrum is unfolded. 
    # The "weights" object stores allele frequencies for each spectrum category
    
    p = seq(0,numChrom)/numChrom
    
    # now get the terms to calculate Tajima's D
    
    a1 = sum(1/seq(1,numChrom))
    a2 = sum(1/seq(1,numChrom)^2)
    b1 = (numChrom+1)/(3*(numChrom-1))
    b2 = 2*(numChrom^2 + numChrom + 3)/(9*numChrom * (numChrom-1))
    c1 = b1 - 1/a1
    c2 = b2 - (numChrom+2)/(a1*numChrom) + a2/a1^2
    e1 = c1/a1
    e2 = c2/(a1^2+a2)
    
    ## Create a dataframe to write the statistics
    out.tab        <- list()
    #names(out.tab) <-c("pop","S","pi","ThetaW")
    out.tab$pop <- pop
    # total number of sites
    
    # Number of variable sites is the sum of the SFS excluding corners
    out.tab$S<-sum(sfs[2:(length(sfs)-1)])
    
    # Calculate Theta as K/an (where an is the harmonic numebr of n-1)
    W<-(out.tab$S/(harmonicNumber))
    out.tab$ThetaW<- W/seq_len
    # calculate pi
    P<-numChrom/(numChrom-1)*2*sum(sfs*p*(1-p))
    out.tab$pi <- P/seq_len
    
    out.tab$TajimaD<-(P - W)/sqrt(e1*out.tab$S+e2*out.tab$S*(out.tab$S-1))
    return(out.tab)
  }))
}

sfs1D <- function(pops, pop,GTs,seq_len){
  n_inds <- length(which(pops==pop))
  
  d <- table(apply(GTs[pops==pop,],2,sum))
  sfs <- d[match(1:(2*n_inds),names(d))]
  
  sfs <- as.vector(c(seq_len-sum(sfs),sfs))
  sfs[is.na(sfs)] <- 0
  return(sfs)
}

f.dist <- function(i, j){
  sqrt((i[1]-j[1])^2 + (i[2]-j[2])^2 )
}


#' finds the origin for a set of Populations
find_origin <- function(psi=tdoa.psi[,5],dist, pct=0.01,xlen=20, ylen=20){
  
  
  dists <- dist[,3]
  ij <- dist[,c(1,2)]
  
  mdlq <- list()
  for(i in 1:xlen){
    mdlq[[i]] <- list()
  }
  
  mdls <- apply(dists,1,function(d){
    l = lm( psi ~ d) 
    l$f.e = .5 * pct/ l$coefficients[2]
    l$rsq = summary(l)$r.squared
    return(l)
  })
  
  lapply(1:nrow(ij),function(r){
    i <- ij[r,1]
    j <- ij[r,2]
    mdlq[[i]][[j]] <<- mdls[[r]]
  })
  
  d0 <- matrix(NA, ncol=xlen, nrow=xlen)
  rsq <- matrix(NA, ncol=xlen, nrow=xlen)
  p <- matrix(NA, ncol=xlen, nrow=xlen)
  
  d0[ij] <- sapply(mdls, function(x)x$f.e)
  rsq[ij] <- sapply(mdls, function(x)x$rsq)
  p[ij] <- sapply(mdls, function(x)summary(x)$coefficients[2,4])
  
  res <- list(d0=d0, rsq=rsq, p, bbox=bbox, xlen=xlen, 
              ylen=ylen)
  class(res) <- 'origin.results'
  
  
  return(res)
  
}


#' gets a bounding box around the populations in samples
get.sample.bbox <- function(samples){
  s <- c('longitude', 'latitude')
  mins <- apply(samples[,s],2,min)
  maxs <- apply(samples[,s],2,max)
  return(cbind(mins, maxs))
}


#' preps data for tdoa analysis
#'
#' makes a data structure of format xi, yi, xj, yj, psi
#' 
#' @param pop.coords coordination data with latitude and longitude
#' @param all.psi matrix with psi values
#' @param region region or set of regions to run analysis on
#' @param countries countries to restrict analysis to
#' @param xlen number of x points
#' @param ylen number of y points
#'
prep.tdoa.data <- function(coords, psi){
  locs <- coords[,c('longitude', 'latitude')]
  n.locs <- nrow(locs)
  
  tdoa.data <- c()
  for(i in 1:n.locs){
    for(j in (i+1):n.locs){
      if( i>=j | j >n.locs) break
      tdoa.data <- rbind( tdoa.data, c(locs[i,], locs[j,], psi[i,j]))
    }
  }
  
  tdoa.data <- matrix(unlist(tdoa.data), ncol=5)
  
  return(tdoa.data)
}



#' Finds origin for a region
#' @param tdoa.data Object of type tdoa.data
#' @param bbox bounding box describing the location where inference should
#'   be done in
#' @param pop.coords coordinations
#' @param pct threshold parameter in model.1d
#' @param xlen, ylen parameters describing the number of points to use
#' @param exclude.ocoean boolean, whether points not on land should be
#'    excluded
#' @param exclude.land boolean, whether points on land should be
#'    excluded
single.origin <- function(tdoa.data, bbox,  pop.coords,
                          pct=0.01,
                          xlen=100, ylen=100, 
                          exclude.ocean=T,
                          exclude.land=F,
                          ...){
  
  
  #define locs for estimate
  s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
  s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
  coords <- expand.grid(s1,s2)
  ij <- expand.grid(1:length(s1), 1:length(s2))
  
  if(exclude.ocean){
    cc <- coords2country(coords)
    to.keep <- !is.na(cc)
  } else {
    to.keep <- rep(T, nrow(coords))
  }
  if(exclude.land){
    cc <- coords2country(coords)
    to.keep <- is.na(cc)
  }
  # init output
  d0 <- matrix(NA, ncol=ylen, nrow=xlen)
  rsq <- matrix(NA, ncol=ylen, nrow=xlen)
  mdlq <- list()
  for(i in 1:xlen){
    mdlq[[i]] <- list()
  }
  
  
  for(r in 1:nrow(coords)){
    i <- ij[r,1]
    j <- ij[r,2]
    x <- coords[r,1]
    y <- coords[r,2]
    
    if(to.keep[r]){
      mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
      d0[i, j] <- mdl$f.e
      rsq[i, j] <- mdl$rsq
      mdlq[[i]][[j]] <- mdl
    }
  }
  res <- list( d0=d0, rsq=rsq, mdlq=mdlq, bbox=bbox, xlen=xlen, 
               ylen=ylen, coords=pop.coords)
  class(res) <- 'origin.results'
  return(res)
}

# tdoa.data <- tdoa.psi
# coords <- pop.coords
single.origin_mod <- function(tdoa.data, bbox,  pop.coords,
                          pct=0.01,
                          xlen=100, ylen=100, 
                          exclude.ocean=F,
                          exclude.land=F,
                          ...){
  
  
  #define locs for estimate
  s1<-seq(bbox[1,1],bbox[1,2], length.out=xlen)
  s2<-seq(bbox[2,1],bbox[2,2], length.out=ylen)
  coords <- expand.grid(s1,s2)
  ij <- expand.grid(1:length(s1), 1:length(s2))
  
  if(exclude.ocean){
    cc <- coords2country(coords)
    to.keep <- !is.na(cc)
  } else {
    to.keep <- rep(T, nrow(coords))
  }
  if(exclude.land){
    cc <- coords2country(coords)
    to.keep <- is.na(cc)
  }
  # init output
  d0 <- matrix(NA, ncol=ylen, nrow=xlen)
  rsq <- matrix(NA, ncol=ylen, nrow=xlen)
  p_val <- matrix(NA, ncol=ylen, nrow=xlen)
  
  # for(i in 1:xlen){
  #   mdlq[[i]] <- list()
  # }
  # 
  r <- 1
  for(r in 1:nrow(coords)){
    i <- ij[r,1]
    j <- ij[r,2]
    x <- coords[r,1]
    y <- coords[r,2]
    
    if(to.keep[r]){
      mdl <- model.1d( xy=c(x,y), data=tdoa.data, pct=pct)
      #summary(mdl)
      d0[i, j] <- mdl$f.e
      rsq[i, j] <- mdl$rsq
      p_val[i, j] <- summary(mdl)$coef[2,4]
      #mdlq[[i]][[j]] <- mdl
    }
  }
  res <- list( d0=d0, rsq=rsq, p_val=p_val, bbox=bbox, xlen=xlen, 
               ylen=ylen, coords=pop.coords)
  class(res) <- 'origin.results'
  return(res)
}


#' gets a bounding box around the populations in samples
get.sample.bbox <- function(samples){
  s <- c('longitude', 'latitude')
  mins <- apply(samples[,s],2,min)
  maxs <- apply(samples[,s],2,max)
  return(cbind(mins, maxs))
}


#' calculates distance from a given point
#' @param xy coordinate of point to evalute function at
#' @param data a 5 col data frame with columns xi, yi, xj, yj, psi, with 
#'   coords from the two sample location and their psi statistic.
#'   best generated unsg prep.tdoa.data
#' @param f.dist the distance function to use. the default 'haversine'
#'   uses the haversine distance. Alternatively, 'euclidean' uses 
#'   Euclidean distance
#' @return an object of type lm describing fit
#' 

#model.1d( xy=c(x,y), data=tdoa.psi, pct=pct) 
model.1d <- function(xy, data, pct=0.01, f.dist="haversine"){
  if (f.dist=="euclidean"){
    f.dist <- function(i, j){
      sqrt((ix -jx)^2 + (iy-jy)^2 )
    }
  }else{ if(f.dist=="haversine"){
    f.dist <- distHaversine
  }}
  
  y = xy[2] 
  x = xy[1]
  ixy = data[,1:2]
  jxy = data[,3:4]
  psi = data[,5]
  
  
  d = f.dist(ixy, c(x,y)) - f.dist(jxy, c(x,y))
  l = lm( psi ~ d ) 
  l$f.e = .5 * pct/ l$coefficients[2]
  l$rsq = summary(l)$r.squared
  
  return (l)
}



#' plots the output of find.origin as a heatmap using .filled.contour
#' @param x an object of type origin.results, as obtained by 
#' @param n.levels the number of color levels
#' @param color.function a function that takes an integer argument and
#'    returns that many colors
#' @param color.negative a single color to be used for negative values
#' @param add.map boolean whether a map should be added
#' @param add.likely.origin boolean, whether origin should be marked with an
#' @param asp aspect ratio, set to 1 to keep aspect ratio with plot
#'   X
#' @export
plot.origin.results <- function(x, n.levels=100, color.function=heat.colors,
                                color.negative='grey',
                                add.map=T,
                                add.samples=T,
                                add.sample.het=T,
                                add.likely.origin=T,
                                asp=1,
                                ...){
  plot.default(NA, xlim=x$bbox[1,], ylim=x$bbox[2,], 
               xlab="", ylab="",
               xaxt="n", yaxt="n",
               xaxs='i', yaxs='i',
               asp=asp, ...)
  
  s1<-seq(x$bbox[1,1],x$bbox[1,2],length.out=x$xlen)
  s2<-seq(x$bbox[2,1],x$bbox[2,2],length.out=x$ylen)
  
  rel <- (x[[1]]>0) * (x[[2]]-min(x[[2]],na.rm=T)) /
    (max(x[[2]],na.rm=T)-min(x[[2]],na.rm=T))+0.001
  
  levels <- c(0,quantile(rel[rel>0.001], 0:n.levels/n.levels, na.rm=T) + 
                1e-6 * 0:n.levels/n.levels)
  cols <- c(color.negative, color.function(n.levels-1))
  
  
  .filled.contour(s1, s2, rel, levels, cols)
  
  
  # rect(x$bbox[1,1], x$bbox[2,1], x$bbox[1,2], x$bbox[2,2], border='black',
  #      lwd=2, col=NULL)
  
  if(add.likely.origin){
    points(summary(x)[,1:2], col='black', pch='x', cex=2)
  }
  if(add.map){
    require(rworldmap)
    m <- getMap("high")
    plot(m, add=T, lwd=1.3)
  }
  
  if(add.sample.het){
    samples <- x$coords
    hets <- (samples$hets - min(samples$hets) )/(
      max(samples$hets) - min(samples$hets))
    points( samples$longitude, samples$latitude,
            pch=16, cex=3, col=grey(hets) )
    points( samples$longitude, samples$latitude,
            pch=1, cex=3, col="black",lwd=2 )
  }
  else if(add.samples){
    points( samples$longitude, samples$latitude,
            pch=16, cex=1, col="black",lwd=1 )
  }
}
#' calculates the psi matrix for a pop object
#' 
#' @param pop population data object from make.pop
#' @param n the sample size which we downsample to
#' @param resampling mode of resampling. Currently, only
#'    hyper is supported
#' @return A n x n matrix of psi values
#' @example examples/example_1.r
#' @export
get.all.psi <- function(pop, n=2,
                        subset=NULL,
                        resampling="hyper"){
  #this function calculates psi for columns i,j, both resampled down to
  # n samples
  
  if(is.null(subset))
    subset <- 1:pop$n
  if( is.logical( subset) )
    subset <- which(subset)
  n.pops <- length(subset)
  mat = matrix( 0, nrow=n.pops, ncol=n.pops )
  for(i in 1:(n.pops-1)){
    for(j in (i+1):n.pops){
      ii <- subset[i]
      jj <- subset[j]
      ni <- unlist(pop$ss[ii,])
      nj <- unlist(pop$ss[jj,])
      fi <- unlist(pop$data[ii,])
      fj <- unlist(pop$data[jj,])
      mat[j,i] <- get.psi( ni, nj, fi, fj, 
                           resampling=resampling, n=n )
      mat[i,j] <- -mat[j,i]
      #print( c(ii, jj))
    }
  }
  
  return(mat)
}

#' the psi statistic calculation
#' the actual calculation of the psi statistic for
#' a single pair of populations. The function requires
#' 4 vectors, all of length equal to the number of snps
#' to be analyzed.
#
#'
#' @param ni the number of sampled haplotypes in population i
#' @param nj the number of sampled haplotypes in population j
#' @param fi the number of derived alleles in population i
#' @param fj the number of derived alleles in population j
#' @param n the number of samples to downsample to
#' @example examples/example_1.r
#' @return psi a matrix of pairwise psi values
get.psi <- function (ni, nj, fi, fj,
                     n=2, resampling="hyper"){
  
  fn <- cbind( fi, ni, fj, nj)
  
  tbl <- table(as.data.frame(fn))
  tbl <- as.data.frame( tbl )
  tbl <- tbl[tbl$Freq > 0,]
  
  tbl$fi <- as.integer(as.character( tbl$fi ))
  tbl$fj <- as.integer(as.character( tbl$fj ))
  tbl$ni <- as.integer(as.character( tbl$ni ))
  tbl$nj <- as.integer(as.character( tbl$nj ))
  
  to.exclude <- tbl$fi == 0 | tbl$fj == 0 | 
    tbl$ni < n | tbl$nj < n
  
  tbl <- tbl[! to.exclude, ]
  
  if(nrow(tbl)==0){ return(NaN)}
  
  
  poly.mat <- matrix(0,nrow=n+1, ncol=n+1)
  poly.mat[2:(n+1),2:(n+1)] <- 1
  poly.mat[n+1,n+1] <- 0
  
  #psi.mat is the contribution to psi for each entry
  psi.mat <- outer(0:n,0:n,FUN=function(x,y)(y-x))
  psi.mat[1,] <- 0
  psi.mat[,1] <- 0
  
  f.contribution <- function(row, b=2){
    a <- 0:b
    f1 <- row[1]
    n1 <- row[2]
    f2 <- row[3]
    n2 <- row[4]
    cnt <- row[5]
    q1 <- choose(b, a) * choose(n1-b, f1-a)/choose(n1,f1)
    q2 <- choose(b, a) * choose(n2-b, f2-a)/choose(n2,f2)
    return( cnt * outer(q2, q1) )
  }
  
  
  resampled.mat <- matrix(rowSums(apply(tbl, 1,
                                        f.contribution)),nrow=n+1)
  
  return( sum(resampled.mat * psi.mat) / sum(resampled.mat * poly.mat) )
}
```

## Historical Demography
```bash
# SMC++
## All
module load BCFtools/1.22-GCC-12.3.0
# module load Miniconda3/23.10.0-1
# source activate /path/to/tools/pyrho_env

outdir=/path/to/smcpp/Centrostephanus_rodgersii
mkdir -p $outdir
cd $outdir

mu_rate=8.6e-9
# mu_rate=5e-9
# mu_rate=1e-9


## input for all populations
## recode "." genotypes to "./."
# cd /path/to/wgs_vcf
# bcftools +fixploidy Cr_uniqID_noRel_filt_HWD_th_noOuts_all.vcf.gz -Oz -o Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz
# tabix Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz

# pop_all=Cr
# bgvcf=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz
# bgvcf_id=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_all.ind


## input for modern samples only
# cd /path/to/wgs_vcf
# bcftools view -S ^Cr_OZhist.txt -Oz -o Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN_fixed.vcf.gz Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz
# tabix Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN_fixed.vcf.gz
# bcftools query -l Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN_fixed.vcf.gz > Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN.ind

pop_all=Cr_allMDN
bgvcf=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN_fixed.vcf.gz
bgvcf_id=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_allMDN.ind

# Inference of the whole population
samples=$(cat $bgvcf_id | paste -sd ,)

## convert vcf to smc for each chrom
# while IFS= read -r chr
# do
#     echo "$chr"
#     singularity run /path/to/tools/smcpp.sif \
#     vcf2smc --cores $SLURM_CPUS_PER_TASK $bgvcf ${pop_all}_${chr}.smc.gz ${chr} ${pop_all}:${samples}
# done < chrom_list.txt

### estimate pop size with all chroms
mkdir -p ${pop_all}_mu${mu_rate}
cd ${pop_all}_mu${mu_rate}
singularity run /path/to/tools/smcpp.sif \
estimate --cores $SLURM_CPUS_PER_TASK --timepoints 1e4 1e8 --knots 40 -o ./ --base ${pop_all}_mu${mu_rate} ${mu_rate} ../${pop_all}_*.smc.gz
mv ${pop_all}_mu${mu_rate}.final.json ../
cd ../
rm -rf ${pop_all}_mu${mu_rate}

### plot inference
singularity run /path/to/tools/smcpp.sif \
plot --cores $SLURM_CPUS_PER_TASK --csv ${pop_all}_mu${mu_rate}_smcpp.pdf ${pop_all}_mu${mu_rate}.final.json

## Meta-population
module load BCFtools/1.22-GCC-12.3.0
# module load Miniconda3/23.10.0-1
# source activate /path/to/tools/pyrho_env

outdir=/path/to/smcpp/Centrostephanus_rodgersii
mkdir -p $outdir
cd $outdir


# mu_rate=8.6e-9
mu_rate=5e-9
# mu_rate=1e-9


## recode "." genotypes to "./."
# cd /path/to/tram/wgs_vcf
# bcftools +fixploidy Cr_uniqID_noRel_filt_HWD_th_noOuts_all.vcf.gz -Oz -o Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz
# tabix Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz

## split by population
# cd /path/to/wgs_vcf
# for pop in NZ OZ RA; do
#   echo ${pop}
#   bcftools view -S Cr_uniqID_noRel_filt_HWD025_recomb_${pop}.ind -Oz -o Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.vcf.gz Cr_uniqID_noRel_filt_HWD_th_noOuts_all_fixed.vcf.gz
#   bcftools query -l Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.vcf.gz > Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.ind
#   tabix Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.vcf.gz
# done

## remove historical samples for OZ
# bcftools view -S ^Cr_OZhist.txt -Oz -o Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn.vcf.gz Cr_uniqID_noRel_filt_HWD_th_noOuts_OZ.vcf.gz
# bcftools query -l Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn.vcf.gz > Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn.ind
# tabix Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn.vcf.gz

# input for metapopulations
species=Cr
# pop=NZ
# pop=OZ
pop=OZmdn
# pop=RA
pop_all="${species}_${pop}"
bgvcf=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.vcf.gz
bgvcf_id=/path/to/wgs_vcf/Cr_uniqID_noRel_filt_HWD_th_noOuts_${pop}.ind


# Inference of the whole population
samples=$(cat $bgvcf_id | paste -sd ,)

## convert vcf to smc for each chrom
# while IFS= read -r chr
# do
#     echo "$chr"
#     singularity run /path/to/tools/smcpp.sif \
#     vcf2smc --cores $SLURM_CPUS_PER_TASK $bgvcf ${pop_all}_${chr}.smc.gz ${chr} ${pop_all}:${samples}
# done < chrom_list.txt

### estimate pop size with all chroms
mkdir -p ${pop_all}_mu${mu_rate}
cd ${pop_all}_mu${mu_rate}
singularity run /path/to/tools/smcpp.sif \
estimate --cores $SLURM_CPUS_PER_TASK --timepoints 1e4 1e8 --knots 40 -o ./ --base ${pop_all}_mu${mu_rate} ${mu_rate} ../${pop_all}_*.smc.gz
mv ${pop_all}_mu${mu_rate}.final.json ../
cd ../
rm -rf ${pop_all}_mu${mu_rate}

### plot inference
singularity run /path/to/tools/smcpp.sif \
plot --cores $SLURM_CPUS_PER_TASK --csv ${pop_all}_mu${mu_rate}_smcpp.pdf ${pop_all}_mu${mu_rate}.final.json

# plot 
library(data.table)
library(ggplot2)

### dataset withuot ouliers
smcpp_files <- list.files("/path/to/smcpp/Centrostephanus_rodgersii",
                          ".csv", full.names = T)

smcpp_res <- do.call(rbind, lapply(smcpp_files[!grepl("mdn", smcpp_files, ignore.case = T)], function(f) {
  s <- fread(f)
  s$mu <- gsub("Cr_.*mu|_smcpp.csv", "", basename(f))
  s$pop <- gsub("Cr_|_mu.+|mu.+", "", basename(f))
  s$pop <- ifelse(s$pop == "", "all", s$pop)
  return(s)
}))

png("/path/to/plot/Centrostephanus_rodgersii/smcpp_combined.png",
    width = 8, height = 5, res = 712, units = "in")
ggplot(smcpp_res) +
  facet_grid(rows = vars(pop)) +
  geom_line(aes(x = x+1, y = y, color = mu)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.border = element_rect(color = "gray90", fill = NA),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("black", "gray50", "gray80")) +
  labs(color = "Mutation rate") +
  xlab("Generations ago") +
  ylab("Ne") +
  scale_x_log10(labels = scales::label_number(), 
                # breaks = c(10^2, 10^4, 10^6, 10^8),
                expand = c(0, 0)) +
  scale_y_log10(labels = scales::label_number(),
                limits = c(1, 200000)) +
  annotation_logticks(color = "gray75", linewidth = 0.2)
dev.off()

### dataset with outliers
smcpp_files <- list.files("/path/to/smcpp/Centrostephanus_rodgersii_outliers",
                          ".csv", full.names = T)

smcpp_res <- do.call(rbind, lapply(smcpp_files[!grepl("mdn", smcpp_files, ignore.case = T)], function(f) {
  s <- fread(f)
  s$mu <- gsub(".*mu|_smcpp.csv", "", basename(f))
  s$pop <- gsub("_mu.+", "", basename(f))
  s$pop <- ifelse(s$pop == "", "all", s$pop)
  return(s)
}))

png("/path/to/plot/Centrostephanus_rodgersii/smcpp_combined_wtOutliers.png",
    width = 8, height = 5, res = 712, units = "in")
ggplot(smcpp_res) +
  facet_grid(rows = vars(pop)) +
  geom_line(aes(x = x+1, y = y, color = mu)) +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        panel.border = element_rect(color = "gray90", fill = NA),
        panel.grid = element_blank()) +
  scale_color_manual(values = c("black", "gray50", "gray80")) +
  labs(color = "Mutation rate") +
  xlab("Generations ago") +
  ylab("Ne") +
  scale_x_log10(labels = scales::label_number(), 
                # breaks = c(10^2, 10^4, 10^6, 10^8),
                expand = c(0, 0)) +
  scale_y_log10(labels = scales::label_number(),
                limits = c(1, max(smcpp_res$y))) +
  annotation_logticks(color = "gray75", linewidth = 0.2)
dev.off()

# StairwayPlot2
cd /nesi/nobackup/ga03714/stairway_plot_v2.2

java -cp /nesi/nobackup/ga03714/stairway_plot_v2.2/stairway_plot_es Stairbuilder CrNZ.blueprint
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrNZ.blueprint.sh
bash CrNZ.blueprint.sh
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrNZ.blueprint.plot.sh
bash CrNZ.blueprint.plot.sh

java -cp /nesi/nobackup/ga03714/stairway_plot_v2.2/stairway_plot_es Stairbuilder CrRA.blueprint
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrRA.blueprint.sh
bash CrRA.blueprint.sh
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrRA.blueprint.plot.sh
bash CrRA.blueprint.plot.sh

java -cp /nesi/nobackup/ga03714/stairway_plot_v2.2/stairway_plot_es Stairbuilder CrAUS.blueprint
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrAUS.blueprint.sh
bash CrAUS.blueprint.sh
sed -i 's/stairway_plot_v2\.1\.2/stairway_plot_v2.2/g' CrAUS.blueprint.plot.sh
bash CrAUS.blueprint.plot.sh

## Example blueprint file
#input setting
popid: CrAUS # id of the population (no white space)
nseq: 170 # number of sequences
L: 570283 # total number of observed nucleic sites, including polymorphic and monomorphic
whether_folded: true # whethr the SFS is folded (true or false) 
SFS: 27937.1778802398 68022.56332536184 91731.16149761213 90005.28937639337 72451.54763791902 53764.31034251985 38958.79614752824 28480.23688155297 21040.93706535772 16204.12495625687 12327.34298102728 9438.074652304836 7286.644536318751 5764.639874919504 4628.740986199497 3735.596051287394 2995.240259856331 2420.748400532977 1975.510112867578 1584.036941643059 1339.001029632048 1097.413483098774 890.995280759069 730.7993727299289 634.1588217485333 557.5541386982705 474.4796608693136 407.6926195437047 355.7216590784132 300.5634078314174 261.9262848790686 230.224442293896 197.5134788207985 172.829899067813 152.7883432713367 127.6517968748399 110.6042910116995 95.59348436490009 86.76116155573499 77.60908427626964 73.71112868152701 67.87960464479801 63.2742465444156 60.63703577827376 53.83916756151967 51.02868081408081 47.31281173240844 44.41703106700109 40.80592740142072 36.54508389212353 31.93342951226092 30.6736109598257 30.93385537151829 29.65908525821475 27.96884531256591 27.33031828770915 23.49915976045667 20.39245753532515 19.02498178182352 18.79755524622849 18.73038811362996 18.22382777500981 18.21331980349738 18.56424752363127 18.28614151093171 16.7737289655312 16.21290750029194 14.32346105804777 13.42175130809196 14.43451464173711 14.64824063887739 13.00235936517615 13.23434649304214 15.18836376207562 16.93624768036435 18.75579390039749 17.65760211376708 17.35438648755281 18.11597540364973 19.60087063184824 20.29213417556884 20.48470961558355 21.54866090776686 22.9501799223302 11.78017545295487 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 # snp frequency spectrum: number of singleton, number of doubleton, etc. (separated by white space)
#smallest_size_of_SFS_bin_used_for_estimation: 1 # default is 1; to ignore singletons, uncomment this line and change this number to 2
largest_size_of_SFS_bin_used_for_estimation: 84 # default is nseq/2 for folded SFS
pct_training: 0.67 # percentage of sites for training
nrand: 42 85 127 166 # number of random break points for each try (separated by white space)
project_dir: CrAUS # project directory
stairway_plot_dir: /nesi/nobackup/ga03714/stairway_plot_v2.1.2/stairway_plot_es # directory to the stairway plot files
ninput: 200 # number of input files to be created for each estimation
theta_upper_bound: 1 # the maximum of theta used in the search algorithm, default is 0.2. Increase this value only if seeing confident intervals converged to the same value as a plateau on the plot.
#dimension_factor: 2000 # this parameter determin the maximum number of iteration of the search algorithm, default is 2000.
#random_seed: 6
#output setting
mu: 8.8e-9 # assumed mutation rate per site per generation
year_per_generation: 1 # assumed generation time (in years)
#plot setting
plot_title: CrAUS_Stairway # title of the plot
xrange: 0, 250 # Time (1k year) range; format: xmin,xmax; "0,0" for default
yrange: 0,0 # Ne (1k individual) range; format: xmin,xmax; "0,0" for default
xspacing: 2 # X axis spacing
yspacing: 2 # Y axis spacing
fontsize: 10 # Font size

# GONE2
script_dir=/path/to/centro/script
work_dir=/path/to/gone/Centrostephanus_rodgersii
mkdir -p $work_dir

## copy the vcf file and script files to the working dir
cd $work_dir

vcf_name=Cr_uniqID_noRel_filt_HWD_th_noOuts
vcf_dir=/path/to/Cr_raw_lcWGS_data/Cr_pop_analysis/thin

## downsample SNPs because GONE can't handle more than 2M SNPs
module load VCFtools/0.1.17-GCC-12.3.0-Perl-5.38.2

THREADS=$SLURM_CPUS_PER_TASK

for pop in all OZ NZ RA;
do
    echo ${pop}
    # vcftools --gzvcf ${vcf_dir}/${vcf_name}_${pop}.vcf.gz \
    # --not-chr Crod1.0_mt \
    # --recode --recode-INFO-all --stdout > ${vcf_name}_${pop}_CHR.vcf
    for rec in 2.5 7.9 12 25;
    do
        echo "rec = ${rec}"
        $script_dir/GONE2/gone2 -g 0 -x -r ${rec} -e -t $THREADS -o ${vcf_name}_${pop}_rec${rec} ${vcf_name}_${pop}_CHR.vcf
    done
done

bcftools view -S ^/path/to/wgs_vcf/Cr_OZhist.txt -Ov -o Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn_CHR.vcf Cr_uniqID_noRel_filt_HWD_th_noOuts_OZ_CHR.vcf

for rec in 2.5 7.9 12 25;
do
    echo "rec = ${rec}"
    $script_dir/GONE2/gone2 -g 0 -x -r ${rec} -e -t $THREADS -o Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn_rec${rec} Cr_uniqID_noRel_filt_HWD_th_noOuts_OZmdn_CHR.vcf
done

cr_gone_files <- list.files("/path/to/gone/Centrostephanus_rodgersii_outliers/", "\\d_GONE2_Ne_mix$",
                            full.names = T)
cr_pop_list <- gsub("Cr_thin_noMT_|_rec.+|rec.+", "", basename(cr_gone_files))
cr_rec_list <- gsub("Cr_thin_noMT_.*rec|_GONE2.+", "", basename(cr_gone_files))

cr_gone_list <- do.call(rbind, lapply(1:length(cr_gone_files), function(i) {
  df <- fread(cr_gone_files[i], skip = 11)
  df$rec <- cr_rec_list[i]
  df$pop <- cr_pop_list[i]
  df
}))

cr_gone_list$rec <- factor(cr_gone_list$rec, c(2.5, 7.9, 12, 25))
png("/path/to/plot/Centrostephanus_rodgersii/gone_g0_x_wtOutliers.png",
    width = 8, height = 5, res = 712, units = "in")
ggplot(cr_gone_list, aes(x = generation, y = Ne_metapop/1000)) +
  facet_grid(rows = vars(rec), cols = vars(pop), labeller = "label_both", scales = "free") +
  geom_line() +
  theme_minimal() +
  theme(strip.text = element_text(size = 14),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        panel.border = element_rect(fill = NA, color = "black")) +
  ylab("Ne (1k individuals)") +
  xlab("Generation ago") +
  scale_y_continuous(limits = c(0, NA))
dev.off()







