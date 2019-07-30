# FastQC
/home/skehoe/scripts/fastqc_report.sh
```
#!/bin/bash 
#SBATCH -J ftest  
#SBATCH -D /home/skehoe/NO_BACKUP/raw_reads/test
#SBATCH -o ftest.%j.out 
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=1 
#SBATCH --mem=50M 
#SBATCH --time=00:20:00 
#SBATCH --mail-type=all 
#SBATCH --mail-user=skehoe@zedat.fu-berlin.de 

RAW_PATH="/home/skehoe/NO_BACKUP/raw_reads/test"
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
#SBATCH -D /home/skehoe/NO_BACKUP/ref
#SBATCH -o quant.%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=09:00:00
#SBATCH --mail-type=all
#SBATCH --mail-user=skehoe@zedat.fu-berlin.de

SALM=/home/perugolate/opt/Salmon-0.8.2_linux_x86_64/bin/salmon

while IFS= read -r i; do
  $SALM quant --index tx_idx \
               --libType A \
               --mates1 /home/skehoe/NO_BACKUP/raw_reads/${i}_R1_001.fastq.gz \
               --mates2 /home/skehoe/NO_BACKUP/raw_reads/${i}_R2_001.fastq.gz \
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
> View(res_AB_up) 1020
> View(res_BA_up) 1206
> View(res_BC_up) 32
> View(res_CB_up) 122
> View(res_AC_up) 614
> View(res_CA_up) 954

* genemap_res_AB_up 1020
* genemap_res_BA_up 1206
* genemap_res_BC_up 32
* genemap_res_CB_up 122
* genemap_res_AC_up 614
* genemap_res_CA_up 954

What about their idx's?
# DAVID
```
DAVID
Functional Annotation
Upload gene list, select identifier, and select list type, submit
Gene list manager: select to limit annoations by one or more species, use all species: Felis catus (14) or unknown (15) You are either not sure which identifier type your list contains, or less than 80% of your list has mapped to your chosen identifier type. Please use the Gene Conversion Tool to determine the identifier type.
    Option 1 (Recommended) continue to submit IDs that DAVID could map
    Option 2 convert gene list
Get Functional Annotation Table
Download file
# res_AB_up - [.txt]
# Import Dataset from Text (base), rename 
> tr_ <- read.delim2("~/Documents/salmon/dirs/tr_.txt")
>   View(tr_)
write.table(tr_, file="GOdavid_AB_up.csv", append = FALSE, sep = "\t", na = "NA", dec = ".", row.names = TRUE, col.names= TRUE)
```
```
data <- read.table("data/Mov10_full_counts.txt", header=T, row.names=1)
meta <- read.table("meta/Mov10_full_meta.txt", header=T, row.names=1)

class(meta)
class(data)

ggplot(XXXX) +
        geom_histogram(aes(x=Sample_1_S1),
        stat="bin", bins=200) +
        x lab = ("Raw expression counts") +
        y lab = ("Number of genes")
        
ggplot(XXXX) +
        geom_histogram(aes(x=Sample_1_S1), 
        stat="bin", bins=200) +
        xlim (-5, 500) +
        xlab ("Raw expression counts")
        ylab ("Number of genes")
        
# for Sample_2_S2 and so on..

PrF_mean_counts <- apply( [,1:3],1,mean)
PrF_variance_counts <- apply( [,1:3],1,var)

PrF_df <- dataframe(PrF_mean_counts, PrF_variance_counts)

# do for PF and SF

ggplot(PrF_df)
        geom_point(aes(x=PrF_mean_counts, v=PrF_variance_counts)) +
        geom_line(aes(x=PrF_mean_counts, y=PrF_mean_counts, color="red")) +
        scale_y_log10() +
        scale_x_log10()
        
# do for PF and SF

vignette("DESeq2")
dds <- DESeqDataSetfromMatrix(countData=data, colData=meta, design=~sampletype)
View(counts(dds))
dds <- estimateSizeFactors(dds)
SizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
write.table(normalized_counts, file="data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)

rld <- rlog(dds, blind=TRUE)
# input matrix of log transformed values
rld <- rlog(dds, blind=TRUE)
rld_mat <- assay(rld)
pca <- prcomp(t(rld_mat))
# create df with metadata + PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df + geom_point(aes(x=PC3, y=PC4, color=sampletype))
# extract the rlog matrix from the object
rld_mat <- assay(rld)
# compute pairwise correlation values
rld_cor <- cor(rld_mat)
head(rld_cor) # check output of cor() make note of the rownames and colnames
pheatmap(rld_cor)
heat.colors <- brewer.pal(6, "Blues")
pheatmap(rld_cor, color=heat.colors, border_color=NA, fontsize=10, fontsize_row=10, height=20)

dds <- DESeq # already created, now run the analysis
dds <- DESeq(dds)
# check size factors
SizeFactors(dds)
# total number of raw counts per sample
colSums(counts(dds), normalized=TRUE))
# now plot dispersion estimates
plotDispEsts(dds)
# now define contrasts, extract results table, and shrink the log2 fold changes

contrast_oe <- c("sampletype", "mov10_overexpression", "control")
res_tableOE_unshrunken <- results(dds, contrast=contrast_oe, alpha=0.05)
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableEO_unshrunken)
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
plotMA(res_tableOE, ylim=c(-2,2))
class(res_tableOE)
mcols(res_table, use.names=TRUE)
res_tableOE %>% data.frame() %>% View()
#define contrasts, extract results table and shrink log2FC
contrast_kd <- c("sampletype", "mov10_knockdown", "control")
res_tableKD <- results(dds, contrast=contrast_kd, alpha=0.05)
res_tableKD <- lfcShrink(dds, contrast=contrast_kd, res=res_tableKD)
# summarize results
summary(res_tableOE)

padj.cutoff <- 0.05
lfc.cutoff <- 0.58

res_tableOE_tb <- res_tableOE %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
sigOE <- res_tableOE_tb %>%
        filter(padj.cutoff & abs(log2FoldChange > lfc.cutoff)
sigOE

# do the same for res_tableKD_tb etc
# create tibbles including row names
mov10_mate <- meta %>% 
        rownames_to_column(var="samplename") %>%
        as_tibble()
normalized_counts <- normalized_counts %>% 
        data.frame() %>%
        rownames_to_column(var="gene") %>%
        as_tibble()
# plot expression for single gene
plotCounts(dds, gene="Mov10", intgroup="sampletype")
# save plotcounts to a dataframe object
d <- plotCounts(dds, gene="Mov10", intgroup="sampletype", returnData=TRUE)
# plotting the MOV10 normalized counts, using samplenames (rownames of d as labels)
ggplot(d, aes(x=sampletype, y=count, color=sampletype)) +
        geom_point(position=position_jitter(w=0.1, h=0)) +
        geom_text_repel(aes(label=rownames(d))) +
        theme_bw() +
        ggtitle("MOV10")
        theme(plot_title=element_text(hjust=0.5)
# order results by padj values
top20_sigOE_genes <- res_tableOE_tb %>%
        arrange(padj) %>%
        pull(gene) %>%
        head(n=20)
# normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
        filter(gene %in% top20_sigOE_genes)
# gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
        gather(colnames(top20_sigOE_norm)[2:9],key="samplename", value="normalized_counts")
# chekc the column header in the "gathered" dataframe
View(gathered_top20_sigOE)
gathered_top20_sigOE <- innerjoin(Mov10_meta, gathered_top20_sigOE)
# plot ggplot2
ggplot(gathered_top20_sigOE) +
        geom_point(aes(x=gene, y=normalized_counts, color=(sampletype)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle(Top 20 significant DEGS") +
        theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        theme(plot.title=element_text(hjust=0.5))
# extract normalized expression for significant genes from the OE and control samples (4:9), set gene column (1) to row names 
norm_ORsig <- normalized_counts[,c(1,4:9)] %>% 
        filter(gene %in% (sigOE$gene) %>%
        data.frame() %>%
        column_to_rownames(var="gene")
# annotate heatmap
annotation <- mov10_meta %>%
        select(samplename, sampletype %>%
        data.frame(row.names="samplename")
# set color palette
heat_colors <- brewer.pal(6, "YlOrRd")
pheatmap(norm_OEsig, 
         color = heat_colors, 
         cluster_rows = T, 
         show_rownames = F,
         annotation = annotation, 
         border_color = NA, 
         fontsize = 10, 
         scale = "row", 
         fontsize_row = 10, 
         height = 20)

# Volcano plot; Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
res_tableOE_tb <- res_tableOE_tb %>% 
                  mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)

ggplot(res_tableOE_tb) +
        geom_point(aes(x = log2FoldChange, y = -log10(padj), colour = threshold_OE)) +
        ggtitle("Mov10 overexpression") +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        #scale_y_continuous(limits = c(0,50)) +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25)))
              
# If we want to know where the top 10 genes (lowest padj) in our DE list are located we could label those dots with the gene name on the Volcano plot using geom_text_repel(). First, we need to order the res_tableOEtibble by padj, and add an additional column to it, to  include on those gene names we want to use to label the plot. 

res_tableOE_tb <- res_tableOE_tb %>% arrange(padj) %>% mutate(genelabels = "")
res_tableOE_tb$genelabels[1:10] <- res_tableOE_tb$gene[1:10]
View(res_tableOE_tb)

ggplot(res_tableOE_tb, aes(x = log2FoldChange, y = -log10(padj))) +
        geom_point(aes(colour = threshold_OE)) +
        geom_text_repel(aes(label = genelabels)) +
        ggtitle("Mov10 overexpression") +
        xlab("log2 fold change") + 
        ylab("-log10 adjusted p-value") +
        theme(legend.position = "none",
              plot.title = element_text(size = rel(1.5), hjust = 0.5),
              axis.title = element_text(size = rel(1.25))) 

# DEGreport can use the DESeq2 results output to make the top20 genes and the volcano plots generated above. 
DEGreport::degPlot(dds = dds, res = res, n = 20, xs = "type", group = "condition") # dds object is output from DESeq2

DEGreport::degVolcano(
    data.frame(res[,c("log2FoldChange","padj")]), # table - 2 columns
    plot_text = data.frame(res[1:10,c("log2FoldChange","padj","id")])) # table to add names
    
# Available in the newer version for R 3.4
DEGreport::degPlotWide(dds = dds, genes = row.names(res)[1:5], group = "condition")

Hypothesis testing: Likelihood ratio test (LRT)
library(DESeq2)
library(DEGreport)

# The full model was specified previously with the `design = ~ sampletype`:
# dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ sampletype)

# Likelihood ratio test
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)

# Extract results
res_LRT <- results(dds_lrt)

# Subset the LRT results to return genes with padj < 0.05
sig_res_LRT <- res_LRT %>%
               data.frame() %>%
               rownames_to_column(var="gene") %>% 
               as_tibble() %>% 
               filter(padj < padj.cutoff)
 
# Get sig gene lists
sigLRT_genes <- sig_res_LRT %>% 
                pull(gene)
                
length(sigLRT_genes)

# Compare to numbers we had from Wald test
nrow(sigOE)
nrow(sigKD)

# Subset results for faster cluster finding (for classroom demo purposes)
clustering_sig_genes <- sig_res_LRT %>%
                  arrange(padj) %>%
                  head(n=1000)


# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)

# What type of data structure is the `clusters` output?
class(clusters)

# Let's see what is stored in the `df` component
head(clusters$df)

# Extract the Group 1 genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%
          filter(cluster == 1)
          
 library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(clusterProfiler)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)

## Explore the grch37 table loaded by the annotables library
grch37

## Return the IDs for the gene symbols in the DE results
idx <- grch37$symbol %in% rownames(res_tableOE)

ids <- grch37[idx, ]

## The gene names can map to more than one Ensembl ID (some genes change ID over time), 
## so we need to remove duplicate IDs prior to assessing enriched GO terms
non_duplicates <- which(duplicated(ids$symbol) == FALSE)

ids <- ids[non_duplicates, ] 

## Merge the IDs with the results 
res_ids <- inner_join(res_tableOE_tb, ids, by=c("gene"="symbol"))   

## Create background dataset for hypergeometric testing using all genes tested for significance in the results                 
allOE_genes <- as.character(res_ids$ensgene)

## Extract significant results
sigOE <- filter(res_ids, padj < 0.05)

sigOE_genes <- as.character(sigOE$ensgene)

## Run GO enrichment analysis 
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "ENSEMBL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
                
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)

write.csv(cluster_summary, "results/clusterProfiler_Mov10oe.csv")

## Dotplot 
dotplot(ego, showCategory=50)

## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms
emapplot(ego, showCategory = 50)

## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange

names(OE_foldchanges) <- sigOE$gene

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)
         
## If some of the high fold changes are getting drowned out due to a large range, you could set a maximum fold change value
OE_foldchanges <- ifelse(OE_foldchanges > 2, 2, OE_foldchanges)
OE_foldchanges <- ifelse(OE_foldchanges < -2, -2, OE_foldchanges)

cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange=OE_foldchanges, 
         vertex.label.font=6)

## Subsetting the ego results without overwriting original `ego` variable
ego2 <- ego

ego2@result <- ego@result[c(1,3,4,8,9),]

## Plotting terms of interest
cnetplot(ego2, 
         categorySize="pvalue", 
         foldChange=OE_foldchanges, 
         showCategory = 5, 
         vertex.label.font=6)
         
## Remove any NA values
res_entrez <- filter(res_ids, entrez != "NA")

## Remove any Entrez duplicates
res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]

## Extract the foldchanges
foldchanges <- res_entrez$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- res_entrez$entrez

## Sort fold changes in decreasing order
foldchanges <- sort(foldchanges, decreasing = TRUE)

head(foldchanges)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = foldchanges, # ordered named vector of fold changes (Entrez IDs are the associated names)
              organism = "hsa", # supported organisms listed below
              nPerm = 1000, # default number permutations
              minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
              pvalueCutoff = 0.05, # padj cutoff value
              verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

## Write GSEA results to file
View(gseaKEGG_results)

write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote=F)

## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03040')

## Output images for a single significant KEGG pathway
detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
pathview(gene.data = foldchanges,
              pathway.id = "hsa03040",
              species = "hsa",
              limit = list(gene = 2, # value gives the max/min limit for foldchanges
              cpd = 1))
              
## Output images for all significant KEGG pathways
get_kegg_plots <- function(x) {
   pathview(gene.data = foldchanges, pathway.id = gseaKEGG_results$ID[x], species = "hsa", 
       limit = list(gene = 2, cpd = 1))
}

purrr::map(1:length(gseaKEGG_results$ID), get_kegg_plots)

# GSEA using gene sets associated with BP Gene Ontology terms
gseaGO <- gseGO(geneList = foldchanges, 
              OrgDb = org.Hs.eg.db, 
              ont = 'BP', 
              nPerm = 1000, 
              minGSSize = 20, 
              pvalueCutoff = 0.05,
              verbose = FALSE) 

gseaGO_results <- gseaGO@result

gseaplot(gseaGO, geneSetID = 'GO:0007423')

biocLite("GSEABase")
library(GSEABase)

# Load in GMT file of gene sets (we downloaded from the Broad Institute [website](http://software.broadinstitute.org/gsea/msigdb/collections.jsp) for MSigDB)

c2 <- read.gmt("/data/c2.cp.v6.0.entrez.gmt.txt")

msig <- GSEA(foldchanges, TERM2GENE=c2, verbose=FALSE)

msig_df <- data.frame(msig)

# Set-up

source("http://bioconductor.org/biocLite.R") 
biocLite("SPIA")
library(SPIA)

## Significant genes is a vector of fold changes where the names are ENTREZ gene IDs. The background set is a vector of all the genes represented on the platform.

background_entrez <- res_entrez$entrez

sig_res_entrez <- res_entrez[which(res_entrez$padj < 0.05), ]

sig_entrez <- sig_res_entrez$log2FoldChange

names(sig_entrez) <- sig_res_entrez$entrez

head(sig_entrez)

spia_result <- spia(de=sig_entrez, all=background_entrez, organism="hsa")

head(spia_result, n=20)

plotP(spia_result, threshold=0.05)

## Look at pathway 03013 and view kegglink
subset(spia_result, ID == "03013")
```
