# FastQC
/home/sxxxxx/scripts/fastqc_report.sh
```
#!/bin/bash 
#SBATCH -J ftest  
#SBATCH -D /home/sxxxxx/NO_BACKUP/raw_reads/test
#SBATCH -o ftest.%j.out 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem=50M 
#SBATCH --time=00:20:00 
#SBATCH --mail-type=all 
#SBATCH --mail-user=

RAW_PATH="/home/sxxxxx/NO_BACKUP/raw_reads/test"
cd $RAW_PATH

for fastq_file in $RAW_PATH; do \
        fastqc $fastq_file \
        ;done 

 rm *.zip 
 mkdir fastqc_reports
 mv *.html fastqc_reports/`
 ```
 # get transcripts and ncRNA
  ```
# get the files
wget ftp://ftp.ensembl.org/pub/release-96/fasta/felis_catus/cds/Felis_catus.Felis_catus_9.0.cds.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-96/fasta/felis_catus/ncrna/Felis_catus.Felis_catus_9.0.ncrna.fa.gz
# combine them
zcat Felis_catus.Felis_catus_9.0.cds.all.fa.gz Felis_catus.Felis_catus_9.0.ncrna.fa.gz > tx.fa 
# make a file to associate transcripts to genes
grep "^>" tx.fa | cut -f1,4 -d" " | sed -e 's/>//g' -e 's/ /,/g' -e 's/gene://g' >> tx2gene.csv
 ```
# index transcripts
```
 salmon index --transcripts tx.fa \
             --index tx_idx
```
# quantify samples
```
$ cat samples.txt
Sample_1_S1
Sample_2_S2
Sample_3_S3
Sample_4_S4
Sample_5_S5
Sample_6_S6
Sample_7_S7
Sample_8_S8
Sample_9_S9
```
```
#!/bin/bash
#SBATCH -J quant
#SBATCH -D /home/sxxxxx/NO_BACKUP/ref
#SBATCH -o quant.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=09:00:00
#SBATCH --mail-type=all
#SBATCH --mail-

SALM=/home/xxxxxx/opt/Salmon-0.8.2_linux_x86_64/bin/salmon

while IFS= read -r i; do
  $SALM quant --index tx_idx \
               --libType A \
               --mates1 /home/sxxxxx/NO_BACKUP/raw_reads/${i}_R1_001.fastq.gz \
               --mates2 /home/sxxxxx/NO_BACKUP/raw_reads/${i}_R2_001.fastq.gz \
               --threads 1 \
               --output $i
done < "samples.txt"
```
```
# Install R
sudo apt update
sudo apt install gdebi libxml2-dev libssl-dev libcurl4-openssl-dev libopenblas-dev r-base r-base-dev

# Install RStudio
cd ~/Downloads
wget https://download1.rstudio.org/rstudio-xenial-1.1.447-amd64.deb
sudo gdebi rstudio-xenial-1.1.447-amd64.deb
printf '\nexport QT_STYLE_OVERRIDE=gtk\n' | sudo tee -a ~/.profile

# Install common packages
R --vanilla << EOF
install.packages(c("tidyverse","data.table","dtplyr","devtools","roxygen2","bit64","readr"), repos = "https://cran.rstudio.com/")
q()
EOF

# Install TDD packages
install.packages("testthis")

# Export to HTML/Excel
R --vanilla << EOF
install.packages(c("htmlTable","openxlsx"), repos = "https://cran.rstudio.com/")
q()
EOF

# Blog tools
R --vanilla << EOF
install.packages(c("knitr","rmarkdown"), repos='http://cran.us.r-project.org')
q()
EOF
sudo apt install python-pip
sudo apt install python3-pip
sudo -H pip install markdown rpy2==2.7.1 pelican==3.7.1
sudo -H pip3 install markdown rpy2==2.9.3 pelican==3.7.1 

# PDF extraction tools
sudo apt install libpoppler-cpp-dev default-jre default-jdk r-cran-rjava
sudo R CMD javareconf
R --vanilla << EOF
library(devtools)
install.packages("pdftools", repos = "https://cran.rstudio.com/")
install_github("ropensci/tabulizer")
q()
EOF

# TTF/OTF fonts usage
sudo apt install libfreetype6-dev
R --vanilla << EOF
install.packages("showtext", repos = "https://cran.rstudio.com/")
q()
EOF

# Cairo for graphic devices
sudo apt install libgtk2.0-dev libxt-dev libcairo2-dev
R --vanilla << EOF
install.packages("Cairo", repos = "https://cran.rstudio.com/")
q()
EOF
```
# load and attach libraries
```
library(readr)
library(tximport)
library(DESeq2)
library(IHW)
library(tidyr)
library(dplyr)
library(cowplot)
library(json)
```
# import transcript counts
```
cat sample_data.csv
sample,type
Sample_1_S1,A
Sample_2_S2,A
Sample_3_S3,A
Sample_4_S4,B
Sample_5_S5,B
Sample_6_S6,B
Sample_7_S7,C
Sample_8_S8,C
Sample_9_S9,C
```
```
# get the names of the directories
dirs <- list.dirs(recursive = FALSE, full.names = FALSE)
# keep the ones that begin with "Sample" (in case there are hidden directories e.g. .git)
dirs <- dirs[grepl("^Sample", dirs)]
# construct paths to salmon output files
files <- file.path(dirs, "quant.sf")
# name the files
names(files) <- dirs
# check that the constructed paths point to files
all(file.exists(files))
```
```
> setwd("/home/shauna/Documents/salmon/dirs")
> dirs <- list.dirs(recursive = FALSE, full.names = FALSE)
> dirs <- dirs[grepl("^Sample", dirs)]
> files <- file.path(dirs, "quant.sf")
> names(files) <- dirs
> all(file.exists(files))
[1] TRUE
> tx2gene <- read_csv("tx2gene.csv")
Error in read_csv("tx2gene.csv") : could not find function "read_csv"
> tx2gene <- read.csv("tx2gene.csv")
> txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
Error in tximport(files, type = "salmon", tx2gene = tx2gene) : 
  could not find function "tximport"
> library(tximport)
> txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
reading in files with read_tsv
1 Error in readInfRepFish(x, type) : 
  importing inferential replicates for Salmon or Sailfish requires package `rjson`.
  to skip this step, set dropInfReps=TRUE
> install.packages("rjson")
> txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
reading in files with read_tsv
1 2 3 4 5 6 7 8 9 
transcripts missing from tx2gene: 1
summarizing abundance
summarizing counts
summarizing length
```
```
# read in data.frame linking transcript id to gene id
tx2gene <- read_csv("tx2gene.csv")
# import the salmon output files
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
```
# import and create sample info
```
# make a data frame containing a single column named "sample"
sample_df <- data.frame(sample = colnames(txi$counts), stringsAsFactors = FALSE)
# here we derive the metadata from file names
sampleData <- read_csv("sample_data.csv")
# join the two in a new object
sampleTable <- dplyr::inner_join(sample_df, sampleData)
# add the sample names as rownames
rownames(sampleTable) <- sampleTable$sample
# double check that the names are in the same order as the counts matrix
all(rownames(sampleTable) == colnames(txi$counts))
```
```
> sample_df <- data.frame(sample = colnames(txi$counts), stringsAsFactors = FALSE)
> sampleData <- read_csv("sample_data.csv")
Error in read_csv("sample_data.csv") : could not find function "read_csv"
> sampleData <- read.csv("sample_data.csv")
> sampleTable <- dplyr::inner_join(sample_df, sampleData)
Joining, by = "sample"
Warning message:
Column `sample` joining character vector and factor, coercing into character vector 
> rownames(sampleTable) <- sampleTable$sample
> all(rownames(sampleTable) == colnames(txi$counts))
[1] TRUE
```
# create DESeqDataSet
```
dds <- DESeqDataSetFromTximport(txi, sampleTable, ~type)
# run DESeq2 on both versions
dds <- DESeq(dds, fitType = "local")
dds_p <- DESeq(dds, fitType = "parametric")
# check the dispersion plots
png("disp_plot.png", width = 2*480)
par(mfrow = c(1,2))
plotDispEsts(dds, main = "local fit")
plotDispEsts(dds_p, main ="parametric fit")
dev.off()
```
```
> library(DESeq2)
> dds <- DESeqDataSetFromTximport(txi, sampleTable, ~type)
using counts and average transcript lengths from tximport
> dds <- DESeq(dds, fitType = "local")
estimating size factors
using 'avgTxLength' from assays(dds), correcting for library size
estimating dispersions
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
> dds_p <- DESeq(dds, fitType = "parametric")
using pre-existing normalization factors
estimating dispersions
found already estimated dispersions, replacing these
gene-wise dispersion estimates
mean-dispersion relationship
final dispersion estimates
fitting model and testing
```
```
> colSums(counts(dds))
Sample_1_S1 Sample_2_S2 Sample_3_S3 Sample_4_S4 Sample_5_S5 Sample_6_S6 Sample_7_S7 Sample_8_S8 Sample_9_S9 
   21718486    29688127     1055911    49308877    27360455    32492627    34705628     7860648    13177842
```
# contrasts
alpha 0.05 padj 0.05 absLFC >=1
```
# this is all the results (significant and non-significant)
res_AB <- results(dds, alpha = 0.05, filterFun = ihw, contrast = c("type", "A", "B"))
# this is just the genes with at least a 2-fold difference in expression
res_AB_sig <- subset(res_AB, padj < 0.05 & abs(log2FoldChange) >= 1)
# abs() is used above to ignore the sign on the fold change value
# i.e. it returns both up- and down-regulated genes in res_AB_sig
# instead you might want the up- and down-regulated genes separately
res_AB_up <- subset(res_AB, padj < 0.05 & log2FoldChange >= 1)
res_AB_dn <- subset(res_AB, padj < 0.05 & log2FoldChange <= -1)
# double check that the numbers add up
(nrow(res_AB_up) + nrow(res_AB_dn)) == nrow(res_AB_sig)
```
> View(res_AB_up) 1020
> View(res_BA_up) 1206
> View(res_BC_up) 32
> View(res_CB_up) 122
> View(res_AC_up) 614
> View(res_CA_up) 954
```
> library(IHW)
> res_AB <- results(dds, alpha = 0.05, filterFun = ihw, contrast = c("type", "A", "B"))
> res_AB_sig <- subset(res_AB, padj < 0.05 & abs(log2FoldChange) >= 1)
> res_AB_up <- subset(res_AB, padj < 0.05 & log2FoldChange >= 1)
> res_AB_dn <- subset(res_AB, padj < 0.05 & log2FoldChange <= -1)
> (nrow(res_AB_up) + nrow(res_AB_dn)) == nrow(res_AB_sig)
[1] TRUE
```
# MA plots
```
# res_AB all (significant and non-significant) 
png("ma_plot_ab.png")
plotMA(res_AB, ylim = c(-10,10), main = "primordial versus primary (all)")
dev.off()
# res_AB_sig genes with at least a 2-fold difference in expression
#res_AB_up up-regulated genes
# res_AB_dn down-regulated genes
```
# PCA
```
rld <- rlog(dds)
pca1 <- plotPCA(rld, intgroup = "type") + ggtitle("500 genes")
pca2 <- plotPCA(rld, intgroup = "type", ntop = nrow(txi$counts)) + ggtitle("All genes")
png("pca_plot.png", width = 2*480)
plot_grid(pca1, pca2)
dev.off()
```
```
#get the transformed counts matrix out as a data frame
rldf <- assay(rld) %>% as.data.frame
#transpose the counts matrix and perform PCA
#scale = FALSE because data are scaled by rlog transform already
pca <- prcomp(t(rldf), scale = FALSE)
#get the PC data out
pcax <- as.data.frame(pca$x)
#add a column for type
pcax$type <- sampleTable$type
#get the percent variance explained
percentVar <- round(100 * (pca$sdev)^2 / sum(pca$sdev^2))
pca12 <- ggplot(pcax, aes(x = PC1, y = PC2, colour = type)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_point(size = 4) + ggtitle("PC1 vs PC2")
pca13 <- ggplot(pcax, aes(x = PC1, y = PC3, colour = type)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC3: ", percentVar[3], "% variance")) +
  geom_point(size = 4) + ggtitle("PC1 vs PC3")
png("pca_plot2.png", width = 2*480)
plot_grid(pca12, pca13)
dev.off()
```
# extracting transformed values
```
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 3)
```
```
BiocInstaller::biocLite("vsn")
library("vsn")
meanSdPlot(assay(ntd))
BiocInstaller::biocLite("hexbin")
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(rld))
meanSdPlot(assay(vsd))
```
```
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))

meanSdPlot(assay(vsd))

meanSdPlot(assay(rld))
```
# Heatmap
```
sampleDists <- dist(t(assay(vsd)))
> library("RColorBrewer")
> sampleDistMatrix <- as.matrix(sampleDists)
> rownames(sampleDistMatrix) <- paste(vsd$type, sep="-")
> colnames(sampleDistMatrix) <- NULL
> colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
> pheatmap(sampleDistMatrix,
+          clustering_distance_rows=sampleDists,
+          clustering_distance_cols=sampleDists,
+          col=colors)
```
# Euclidian distances
```
rld_dds <- rlog(dds, fitType = "local")
sampleDists_salm <- dist(t(assay(rld_dds))) #apply the dist function to the transpose of the transformed count matrix to get sample-to-sample distances
sampleDistMatrix_salm <- as.matrix(sampleDists_salm)
pheatmap(sampleDistMatrix_salm, clustering_distance_rows = sampleDists_salm, clustering_distance_cols = sampleDists_salm, main = "Euclidian distances between the samples")
```
# pheatmaps
```
library(DESeq2)
select_salm <- order(rowMeans(counts(dds, normalized=T)),decreasing = T)[1:50]
nt_salm <- normTransform(dds)
log2.norm.counts_salm <- assay(nt_salm)[select_salm,]
install.packages("pheatmap")
pheatmap(log2.norm.counts_salm, cluster_rows = F, show_rownames = T, cluster_cols = T, legend = T, main = "Pheatmap 50 most highly expressed genes")
Error in pheatmap(log2.norm.counts_salm, cluster_rows = F, show_rownames = T,  : 
  could not find function "pheatmap"
library(pheatmap)
pheatmap(log2.norm.counts_salm, cluster_rows = F, show_rownames = T, cluster_cols = T, legend = T, main = "Pheatmap 50 most highly expressed genes")

(log2.norm.counts_salm)
```
# adding gene names
```
res_AC_up$ensembl <- sapply(strsplit(rownames(res_AC_up), split="\\+"), "[", 1)
library("biomaRt")
ensembl = useMart("ensembl", dataset = "fcatus_gene_ensembl")
genemap_res_AC_up <- getBM(attributes = c("ensembl_gene_id", "ensembl_gene_id_version", "entrezgene", "hgnc_symbol"),
                  filters = "ensembl_gene_id_version",
                  values = res_AC_up$ensembl,
                  mart = ensembl)
idx_res_AC_up <- match(res_AC_up$ensembl, genemap_res_AC_up$ensembl_gene_id_version)
```
# plotting genes
```
library(DESeq2)
library(dplyr)
library(tidyr)
library(MASS)
library(visreg)
library(cowplot)
# get the normalized counts out
ddc <- counts(dds, normalized = TRUE) %>% as.data.frame
# carry over the gene names
ddc$gene_id <- rownames(ddc)
# reformat from wide to long
ddc_long <- gather(ddc, key = library, value = count, -gene_id)
# add the sample data
ddc_long_meta <- inner_join(ddc_long, data.frame(dds@colData, library = dds@colData[,1]))
# This is a bit backwards because it involves taking the normalized counts back out, fitting a GLM, 
# and then extracting the coefficients and intervals for plotting. 
# DESeq has already fitted GLMs for each gene but the data aren't in a very useful state.
# plot a given gene
glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000028290.3"))) %>%
  visreg
bmp15_fit <- glm.nb(count ~ type, data = subset(ddc_long_meta, gene_id %in% c("ENSFCAG00000028290.3"))) %>%
  visreg(plot = FALSE) %$% fit %>% dplyr::select(-count)
# visreg returns logged coefficients/bounds so compute the exponentials
bmp15_fit$visregFit <- exp(bmp15_fit$visregFit)
bmp15_fit$visregUpr <- exp(bmp15_fit$visregUpr)
bmp15_fit$visregLwr <- exp(bmp15_fit$visregLwr)
# 
ggplot(bmp15_fit, aes(x = type, y = visregFit)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), limits = c(0.1, 100000)) +
  geom_pointrange(aes(ymin = visregLwr, ymax = visregUpr), size = 1) +
  xlab("Follicle type") + ylab("Gene expression (normalized count)")
# With n=3 it's probably worth putting the actual data points on too.
# DESeq2 has already fitted glms. So it would be possible to use those values instead. 
# The coefficients it provides aren't in a very useful format:
coef(dds, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3")
# So if you want to use them directly from DESeq, I guess the easiest thing to do would be to remove the intercept. i.e. run it again without an intercept:
dds2 <- dds
dds2@design <- ~ type + 0
dds2 <- DESeq(dds2)
# This would give more useful coefficients:
coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3")
# but they are log2, so:
2^coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t
# These are almost identical to the coefficients from glm.nb:
bmp15_fit
# The intervals could will be a bit different for some genes (e.g. because of the shrinkage of the dispersion parameter in DESeq). 
# For "ENSFCAG00000028290.3", they are basically the same.
ce <- coef(dds2, SE = FALSE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t %>% as.data.frame
cee <- coef(dds2, SE = TRUE) %>% subset(., rownames(.) == "ENSFCAG00000028290.3") %>% t %>% as.data.frame
bmp15_fd <- data.frame(type = rownames(ce), coef = ce[,1], se = cee[,1])
bmp15_fd$lower <- bmp15_fd$coef - (2*bmp15_fd$se)
bmp15_fd$upper <- bmp15_fd$coef + (2*bmp15_fd$se)
# From DESeq:
bmp15_fd %>% dplyr::select(-se) %>% mutate_at(c("coef", "lower", "upper"), function(x) (2^x))
# From basic glm:
bmp15_fit
```

	

