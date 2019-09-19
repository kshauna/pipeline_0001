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
# Summary 
* DESeq2 output
       * 0.05_1 ENSFCAGxx baseMean log2FC lfcSE stat pvalue padj weight
       * res_AB_sigord - up_ord and dn_ord + AC and BC, respectively

I have the following outputs:
* **bioMart genemap output:** No. ensembl_gene_id ensembl_gene_id_version entrez_gene_id hgnc_symbol
       * results for genemap_AB - up and dn + AC and BC, respectively
       * **to do** write.table for NA values sort() arrange()
       * **to do** write.table take out NA values, take out No. and alphabetically sort by hgnc sort() arrange() 

* **bioMart gores output:** No. emsembl_gene_id ensembl_gene_id_version go_id name_1006 def_1006
       * results for gores_AB - up and dn + AC and BC, respectively
       * **to do** write.table for NA values sort() arrange()
       * **to do** write.table take out NA values, take out No. and alpha-numerically sort by go_id sort() arrange()

* **DESeq2 results** for res_AB - up and dn + AC and BC, respectively + log2FC coloured
       * **to do** write.table numerically sort by log2FC sort() arrange() highest-lowest

* **DAVID genelist output:** entrez_gene_ID Name Species
       * results for genelist_AB - up and dn + AC and BC, respectively
       * **to do** .txt to .csv
       * **to do** write.table take out NA values, take out Species and alpha sort by Name sort() arrange()

* **DAVID func. annot. table output:** ID Gene Name Species COG ont GOBP GOCC GO MF INTERPRO KEGG PATH PIR SUPERFAM SMART UP KEY UP SEQ FEATS..... 
      * results for AB dn only **should I continue with DAVID?**
      * **to do** the rest
      * **to do** write.table take out NA values, Species, COG ont INTERPRO KEGG PATH PIR SUPERFAM SMART UP KEY UP SEQ FEATS....., create separate ID Gene Name GOBP GOCC GO MF, and sort alphabetically

* **bioMart querychromo output:** No. ensembl_gene_id ensembl_gene_id_version hgnc chromosome_name stat position end position band
      * results for querychromo_AB - up and dn + AC and BC, respectively
       
What other outputs do I need? 
* **bioMart pvalue output:** hgnc p_value
      * I need this as input for topGO
      * needs to be padjs from res_AB - up and dn + AC and BC, respectively
      * **left col** hgnc (or gene ID) **right col** p_value (padj)
      * results for pvalGO_AB - up and dn + AC and BC, respectively
      * **output** ggplot(goEnrichment..)

# topGO
```
write.table(res_BC,"res_BC.csv",sep=",",quote=FALSE,row.names=TRUE, col.names = NA)
write.table(res_BC_sig,"res_BC_sig.csv",sep=",",quote=FALSE,row.names=TRUE, col.names = NA)

# Gene universe file
BC_univ_exp_data <- read.table('res_BC.csv', header=TRUE, sep=',')
BC_univ_genes <- as.character(BC_univ_exp_data[,1])

# Read in genes of interest
BC_sig_candi_list <- read.table('res_BC_sig.csv', header=TRUE, sep=',')
BC_sig_candi_list <- as.character(BC_sig_candi_list[,1])

# Step 2: Create GO annotation
# create GO db for genes to be used using biomaRt
db <- useMart('ENSEMBL_MART_ENSEMBL', dataset='fcatus_gene_ensembl', host="www.ensembl.org")

BCuniv_go_ids <- getBM(attributes=c('go_id', 'ensembl_gene_id_version', 'namespace_1003'), 
                       filters='ensembl_gene_id_version', 
                       values=BC_univ_genes, mart=db)

# build the gene 2 GO annotation list (needed to create topGO object)
BCuniv_gene2GO <- unstack(BCuniv_go_ids[,c(1,2)])

# make named factor showing which genes are of interest
BC_geneList <- factor(as.integer(BC_univ_genes%in%BC_sig_candi_list))
names(BC_geneList) <- BC_univ_genes

# Step 3: Make topGO data object
BC_GOdata <- new('topGOdata', ontology='BP', 
                 allGenes = BC_geneList, 
                 annot = annFUN.gene2GO, 
                 gene2GO = BCuniv_gene2GO)

# define test using the weight01 algorithm (default) with fisher
weight_fish_resBC <- runTest(BC_GOdata, algorithm='weight01', statistic='fisher') 

# generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
BCallGO <- usedGO(BC_GOdata)
all_resBC <- GenTable(BC_GOdata, 
                      weightFisherBC=weight_fish_resBC, 
                      orderBy='weightFisher', 
                      topNodes=length(BCallGO))

# Step 5: Correcting for multiple testing
#performing BH correction on our p values
p.adj <- round(p.adjust(all_resBC$weightFisherBC,method="BH"), digits = 4)

# create the file with all the statistics from GO analysis
all_resBC_final <- cbind(all_resBC,p.adj)
all_resBC_final <- all_resBC_final[order(all_resBC_final$p.adj),]

#get list of significant GO before multiple testing correction
resBC.table.p <- all_resBC_final[which(all_resBC_final$weightFisherBC<=0.001),]

#get list of significant GO after multiple testing correction
resBC.table.bh <- all_resBC_final[which(all_resBC_final$p.adj<=0.05),]

#save first top 50 ontolgies sorted by adjusted pvalues
write.table(all_resBC_final[1:50,], "summary_BCtopGO_analysis.csv",sep=",", quote=FALSE, row.names=FALSE)

# Plot the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level

pdf(file='topGOPlotBC_fullnames.pdf', 
    height=12, width=12, 
    paper='special', 
    pointsize=18)

showSigOfNodes(BC_GOdata, score(weight_fish_resBC), useInfo = "none", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
dev.off()

# Step 6: Get all the genes in your significant GO TERMS

myBCterms <- resBC.table.p$GO.ID # change it to resBC.table.bh$GO.ID if working with BH corrected values
myBCgenes <- genesInTerm(BC_GOdata, myBCterms)

var=c()
for (i in 1:length(myBCterms))
{
  myBCterm=myBCterms[i]
  myBCgenesforterm= myBCgenes[myBCterm][[1]]
  myBCgenesforterm=paste(myBCgenesforterm, collapse=',')
  var[i]=paste("GOTerm",myBCterm,"genes-", myBCgenesforterm)
}

write.table(var,"genetoGOBCmapping.txt",sep="\t",quote=F)

# print the package versions used ---#
sessionInfo()
```
# REVIGO
Scatterplot & Table
```
  # Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
  # terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800
  
  
  # install.packages( "ggplot2" );
  library( ggplot2 );

  # install.packages( "scales" );
  library( scales );
  
  revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
  revigo.data <- rbind(c("GO:0007160","cell-matrix adhesion", 0.051, 3.906,-4.903, 3.817,-3.1024,0.911,0.000),
  c("GO:0007229","integrin-mediated signaling pathway", 0.056, 2.806, 6.449, 3.860,-5.3872,0.736,0.000),
  c("GO:0048251","elastic fiber assembly", 0.001,-2.279,-4.957, 2.127,-4.4437,0.779,0.043),
  c("GO:0043046","DNA methylation involved in gamete generation", 0.004,-5.420, 2.138, 2.663,-3.3098,0.667,0.046),
  c("GO:0046839","phospholipid dephosphorylation", 0.052, 1.157,-2.521, 3.825,-2.2240,0.834,0.098),
  c("GO:0034587","piRNA metabolic process", 0.004,-5.450,-2.828, 2.731,-2.4157,0.914,0.119),
  c("GO:0019915","lipid storage", 0.032, 6.014, 2.701, 3.611,-2.4342,0.788,0.162),
  c("GO:0055114","oxidation-reduction process",15.060,-0.646,-3.693, 6.286,-3.0223,0.854,0.166),
  c("GO:0010951","negative regulation of endopeptidase activity", 0.157, 1.507, 2.617, 4.304,-4.7696,0.748,0.194),
  c("GO:1990535","neuron projection maintenance", 0.000,-1.176,-5.343, 1.301,-2.3354,0.821,0.196),
  c("GO:0019800","peptide cross-linking via chondroitin 4-sulfate glycosaminoglycan", 0.001,-5.736,-3.327, 2.260,-2.4413,0.831,0.196),
  c("GO:0030335","positive regulation of cell migration", 0.076, 4.275, 2.717, 3.988,-2.5850,0.625,0.199),
  c("GO:0018158","protein oxidation", 0.002,-6.333,-2.047, 2.398,-2.1798,0.829,0.215),
  c("GO:0071230","cellular response to amino acid stimulus", 0.019, 1.030, 8.797, 3.390,-3.7212,0.848,0.265),
  c("GO:0032489","regulation of Cdc42 protein signal transduction", 0.001, 1.965, 6.841, 2.155,-2.4413,0.769,0.268),
  c("GO:0051382","kinetochore assembly", 0.018,-2.523,-4.318, 3.354,-3.0915,0.818,0.358),
  c("GO:0035987","endodermal cell differentiation", 0.011,-4.873, 4.065, 3.156,-2.9586,0.646,0.378),
  c("GO:0043932","ossification involved in bone remodeling", 0.001,-5.997, 3.923, 2.064,-2.3354,0.746,0.379),
  c("GO:0051592","response to calcium ion", 0.018, 0.376, 8.954, 3.357,-2.4449,0.885,0.393),
  c("GO:0034375","high-density lipoprotein particle remodeling", 0.003,-3.103, 2.456, 2.568,-2.2211,0.614,0.402),
  c("GO:1900272","negative regulation of long-term synaptic potentiation", 0.000, 3.520, 4.563, 1.462,-2.1791,0.705,0.403),
  c("GO:0035469","determination of pancreatic left/right asymmetry", 0.002,-5.413, 4.850, 2.394,-2.3354,0.716,0.410),
  c("GO:0032964","collagen biosynthetic process", 0.005,-5.499, 3.359, 2.789,-2.1824,0.713,0.412),
  c("GO:0060047","heart contraction", 0.045,-5.517, 3.862, 3.761,-2.4498,0.758,0.435),
  c("GO:0035137","hindlimb morphogenesis", 0.008,-5.528, 4.375, 2.997,-2.1785,0.708,0.442),
  c("GO:0030195","negative regulation of blood coagulation", 0.009,-1.612, 5.012, 3.064,-2.4413,0.583,0.487),
  c("GO:0051603","proteolysis involved in cellular protein catabolic process", 0.759, 1.503,-3.773, 4.988,-2.8069,0.883,0.493),
  c("GO:0060113","inner ear receptor cell differentiation", 0.016,-4.791, 4.424, 3.300,-2.3391,0.649,0.496),
  c("GO:2001046","positive regulation of integrin-mediated signaling pathway", 0.001, 3.341, 5.722, 2.270,-3.6021,0.700,0.499),
  c("GO:0038026","reelin-mediated signaling pathway", 0.002, 3.073, 7.107, 2.436,-2.7696,0.767,0.509),
  c("GO:0048739","cardiac muscle fiber development", 0.001,-5.080, 5.016, 2.279,-2.8633,0.666,0.523),
  c("GO:0048549","positive regulation of pinocytosis", 0.001, 3.361, 1.230, 2.134,-2.3363,0.683,0.524),
  c("GO:0051781","positive regulation of cell division", 0.023, 4.410, 3.960, 3.477,-2.4145,0.725,0.525),
  c("GO:1905342","positive regulation of protein localization to kinetochore", 0.001, 5.471, 2.832, 2.173,-2.3354,0.747,0.527),
  c("GO:0007129","synapsis", 0.018,-2.580,-2.290, 3.362,-3.1675,0.707,0.555),
  c("GO:1901379","regulation of potassium ion transmembrane transport", 0.014, 4.968, 1.625, 3.258,-2.4473,0.698,0.609),
  c("GO:0034113","heterotypic cell-cell adhesion", 0.009, 4.136,-4.744, 3.053,-2.3197,0.852,0.613),
  c("GO:0060686","negative regulation of prostatic bud formation", 0.001,-3.170, 3.846, 2.097,-2.3354,0.573,0.616),
  c("GO:0033627","cell adhesion mediated by integrin", 0.013, 4.400,-4.540, 3.231,-2.1752,0.912,0.628),
  c("GO:0043691","reverse cholesterol transport", 0.003, 7.106, 0.687, 2.646,-2.2211,0.806,0.634),
  c("GO:0060397","JAK-STAT cascade involved in growth hormone signaling pathway", 0.001, 2.210, 7.661, 1.934,-2.3354,0.752,0.684),
  c("GO:0072678","T cell migration", 0.007, 6.438,-0.956, 2.961,-2.3382,0.776,0.695));
  
  one.data <- data.frame(revigo.data);
  names(one.data) <- revigo.names;
  one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
  one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
  one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
  one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
  one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
  one.data$frequency <- as.numeric( as.character(one.data$frequency) );
  one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
  one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
  head(one.data);
  
  
  # --------------------------------------------------------------------------
  # Names of the axes, sizes of the numbers and letters, names of the columns,
  # etc. can be changed below
  
  p1 <- ggplot( data = one.data );
  p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
  p1 <- p1 + scale_colour_gradientn( colours = c("blue", "green", "yellow", "red"), limits = c( min(one.data$log10_p_value), 0) );
  p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
  p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
  ex <- one.data [ one.data$dispensability < 0.15, ]; 
  p1 <- p1 + geom_text( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 3 );
  p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
  p1 <- p1 + theme(legend.key = element_blank()) ;
  one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
  one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
  p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
  p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);
  
  
  
  # --------------------------------------------------------------------------
  # Output the plot to screen
  
  p1;
  
  # Uncomment the line below to also save the plot to a file.
  # The file type depends on the extension (default=pdf).
  
  # ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
```
Treemap
```
# A treemap R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

# author: Anton Kratz <anton.kratz@gmail.com>, RIKEN Omics Science Center, Functional Genomics Technology Team, Japan
# created: Fri, Nov 02, 2012  7:25:52 PM
# last change: Fri, Nov 09, 2012  3:20:01 PM

# -----------------------------------------------------------------------------
# If you don't have the treemap package installed, uncomment the following line:
# install.packages( "treemap" );
library(treemap) 								# treemap package by Martijn Tennekes

# Set the working directory if necessary
# setwd("C:/Users/username/workingdir");

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0007160","cell-matrix adhesion",0.051,3.1024,0.911,0.000,"cell-matrix adhesion"),
c("GO:0033627","cell adhesion mediated by integrin",0.013,2.1752,0.912,0.628,"cell-matrix adhesion"),
c("GO:0034113","heterotypic cell-cell adhesion",0.009,2.3197,0.852,0.613,"cell-matrix adhesion"),
c("GO:0007229","integrin-mediated signaling pathway",0.056,5.3872,0.736,0.000,"integrin-mediated signaling pathway"),
c("GO:0071230","cellular response to amino acid stimulus",0.019,3.7212,0.848,0.265,"integrin-mediated signaling pathway"),
c("GO:0072678","T cell migration",0.007,2.3382,0.776,0.695,"integrin-mediated signaling pathway"),
c("GO:0030335","positive regulation of cell migration",0.076,2.5850,0.625,0.199,"integrin-mediated signaling pathway"),
c("GO:1905342","positive regulation of protein localization to kinetochore",0.001,2.3354,0.747,0.527,"integrin-mediated signaling pathway"),
c("GO:0019915","lipid storage",0.032,2.4342,0.788,0.162,"integrin-mediated signaling pathway"),
c("GO:0051592","response to calcium ion",0.018,2.4449,0.885,0.393,"integrin-mediated signaling pathway"),
c("GO:1900272","negative regulation of long-term synaptic potentiation",0.000,2.1791,0.705,0.403,"integrin-mediated signaling pathway"),
c("GO:0048549","positive regulation of pinocytosis",0.001,2.3363,0.683,0.524,"integrin-mediated signaling pathway"),
c("GO:0051603","proteolysis involved in cellular protein catabolic process",0.759,2.8069,0.883,0.493,"integrin-mediated signaling pathway"),
c("GO:0010951","negative regulation of endopeptidase activity",0.157,4.7696,0.748,0.194,"integrin-mediated signaling pathway"),
c("GO:1901379","regulation of potassium ion transmembrane transport",0.014,2.4473,0.698,0.609,"integrin-mediated signaling pathway"),
c("GO:0043691","reverse cholesterol transport",0.003,2.2211,0.806,0.634,"integrin-mediated signaling pathway"),
c("GO:2001046","positive regulation of integrin-mediated signaling pathway",0.001,3.6021,0.700,0.499,"integrin-mediated signaling pathway"),
c("GO:0060397","JAK-STAT cascade involved in growth hormone signaling pathway",0.001,2.3354,0.752,0.684,"integrin-mediated signaling pathway"),
c("GO:0038026","reelin-mediated signaling pathway",0.002,2.7696,0.767,0.509,"integrin-mediated signaling pathway"),
c("GO:0032489","regulation of Cdc42 protein signal transduction",0.001,2.4413,0.769,0.268,"integrin-mediated signaling pathway"),
c("GO:0030195","negative regulation of blood coagulation",0.009,2.4413,0.583,0.487,"integrin-mediated signaling pathway"),
c("GO:0051781","positive regulation of cell division",0.023,2.4145,0.725,0.525,"integrin-mediated signaling pathway"),
c("GO:0048251","elastic fiber assembly",0.001,4.4437,0.779,0.043,"elastic fiber assembly"),
c("GO:1990535","neuron projection maintenance",0.000,2.3354,0.821,0.196,"elastic fiber assembly"),
c("GO:0051382","kinetochore assembly",0.018,3.0915,0.818,0.358,"elastic fiber assembly"),
c("GO:0043046","DNA methylation involved in gamete generation",0.004,3.3098,0.667,0.046,"DNA methylation involved in gamete generation"),
c("GO:0018158","protein oxidation",0.002,2.1798,0.829,0.215,"DNA methylation involved in gamete generation"),
c("GO:0007129","synapsis",0.018,3.1675,0.707,0.555,"DNA methylation involved in gamete generation"),
c("GO:0019800","peptide cross-linking via chondroitin 4-sulfate glycosaminoglycan",0.001,2.4413,0.831,0.196,"DNA methylation involved in gamete generation"),
c("GO:0034375","high-density lipoprotein particle remodeling",0.003,2.2211,0.614,0.402,"DNA methylation involved in gamete generation"),
c("GO:0060686","negative regulation of prostatic bud formation",0.001,2.3354,0.573,0.616,"DNA methylation involved in gamete generation"),
c("GO:0060113","inner ear receptor cell differentiation",0.016,2.3391,0.649,0.496,"DNA methylation involved in gamete generation"),
c("GO:0034587","piRNA metabolic process",0.004,2.4157,0.914,0.119,"DNA methylation involved in gamete generation"),
c("GO:0043932","ossification involved in bone remodeling",0.001,2.3354,0.746,0.379,"DNA methylation involved in gamete generation"),
c("GO:0035469","determination of pancreatic left/right asymmetry",0.002,2.3354,0.716,0.410,"DNA methylation involved in gamete generation"),
c("GO:0035987","endodermal cell differentiation",0.011,2.9586,0.646,0.378,"DNA methylation involved in gamete generation"),
c("GO:0032964","collagen biosynthetic process",0.005,2.1824,0.713,0.412,"DNA methylation involved in gamete generation"),
c("GO:0048739","cardiac muscle fiber development",0.001,2.8633,0.666,0.523,"DNA methylation involved in gamete generation"),
c("GO:0035137","hindlimb morphogenesis",0.008,2.1785,0.708,0.442,"DNA methylation involved in gamete generation"),
c("GO:0060047","heart contraction",0.045,2.4498,0.758,0.435,"DNA methylation involved in gamete generation"),
c("GO:0046839","phospholipid dephosphorylation",0.052,2.2240,0.834,0.098,"phospholipid dephosphorylation"),
c("GO:0055114","oxidation-reduction process",15.060,3.0223,0.854,0.166,"phospholipid dephosphorylation"));

stuff <- data.frame(revigo.data);
names(stuff) <- revigo.names;

stuff$abslog10pvalue <- as.numeric( as.character(stuff$abslog10pvalue) );
stuff$freqInDbPercent <- as.numeric( as.character(stuff$freqInDbPercent) );
stuff$uniqueness <- as.numeric( as.character(stuff$uniqueness) );
stuff$dispensability <- as.numeric( as.character(stuff$dispensability) );

# by default, outputs to a PDF file
pdf( file="revigo_treemap.pdf", width=16, height=9 ) # width and height are in inches

# check the tmPlot command documentation for all possible parameters - there are a lot more
tmPlot(
	stuff,
	index = c("representative","description"),
	vSize = "abslog10pvalue",
	type = "categorical",
	vColor = "representative",
	title = "REVIGO Gene Ontology treemap",
	inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
	lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
	bg.labels = "#CCCCCCAA",     # define background color of group labels
												       # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
	position.legend = "none"
)

dev.off()
```
### 160919
```
> View(res_BA_dn) 1020
> View(res_BA_sig) 2226
> View(res_BA_up) 1206
> View(res_CA_dn) 614
> View(res_CA_sig) 1568
> View(res_CA_up) 954
> View(res_CB_dn) 32 
> View(res_CB_sig) 154 
> View(res_CB_up) 122

> View(res_AB_dn) 1206
> View(res_AB_sig) 2226
> View(res_AB_up) 1020
> View(res_AC_dn) 954
> View(res_AC_sig) 1568
> View(res_AC_up) 614
> View(res_BC_dn) 122
> View(res_BC_sig) 154
> View(res_BC_up) 32
```
### 190919
topGO statistics
* define test using the weight01 algorithm (default) with fisher
	* weight_fish_resBCMF <- runTest(BCMF_GOdata, algorithm='weight01', statistic='fisher') 
	* generate a table of results: using weightFisher & orderBy = weightFisher
	* Correcting for multiple testing performing BH correction on our p values round(p.adjust(all_resBCMF$weightFisherBCMF,method="BH")
	* create the file with all the statistics from GO analysis
	* get list of significant GO before multiple testing correction weightFisher <=0.001
	* get list of significant GO after multiple testing correction p.adj <=0.05
	* write.table save first top 50 ontolgies sorted by adjusted pvalues
	

