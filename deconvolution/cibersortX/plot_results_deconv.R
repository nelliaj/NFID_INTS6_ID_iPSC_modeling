library(tidyverse)
library(stats)
library(RColorBrewer)
library(pheatmap)
library(yaml)

## Pau's reannotation
yml_file <- read_yaml("deconvolution/Data_S1_snapseed_markers_updated.yaml", fileEncoding = "UTF-8")

deconv_results <- read.table("deconvolution/cibersortx/output/snapseed/abbreviation/CIBERSORTx_Adjusted.txt", header=T)

metadata_bulk <- read.table("deconvolution/metadata_with_rnu4-1.txt", header=T)
metadata_bulk <- metadata_bulk[metadata_bulk$tp=="d28",]

matrix_deconv <- deconv_results[,2:11]
rownames(matrix_deconv) <- deconv_results$Mixture
pheatmap(matrix_deconv)


meta_matrix_deconv <- data.frame(Mixture=deconv_results$Mixture)
meta_matrix_deconv$phenotype <- metadata_bulk[match(meta_matrix_deconv$Mixture, metadata_bulk$seqsample),]$phenotype
meta_matrix_deconv$sex <- metadata_bulk[match(meta_matrix_deconv$Mixture, metadata_bulk$seqsample),]$sex
meta_matrix_deconv$astro_percentage <- metadata_bulk[match(meta_matrix_deconv$Mixture, metadata_bulk$seqsample),]$astro_percentage
rownames(meta_matrix_deconv) <- meta_matrix_deconv$Mixture
meta_matrix_deconv$Mixture <- NULL

## Heatmaps of cell type proportion

plot_prop <- pheatmap(
  matrix_deconv,
  annotation_row = meta_matrix_deconv,
  main = "Cell type proportions (cibersortX)",
  fontsize_col = 12,
  fontsize_row=12
)

pdf("deconvolution/plots/heatmap_cTypeproportions.pdf", width = 8, height = 10)
plot_prop
dev.off()


## Stacked barplot of iPSC lines


matrix_deconv_new <- deconv_results[,1:11]
matrix_deconv_new_long <- as.data.frame(matrix_deconv_new %>% pivot_longer(-c("Mixture"), names_to="CellType", values_to="Proportion"))


matrix_deconv_new_collapsed <- deconv_results[,1:11]

en_collapsed <- rowSums(matrix_deconv_new_collapsed[,grepl("EN", colnames(matrix_deconv_new_collapsed))])
matrix_deconv_new_collapsed <- matrix_deconv_new_collapsed[,!grepl("EN", colnames(matrix_deconv_new_collapsed))]
matrix_deconv_new_collapsed$EN <- en_collapsed

npc_collapsed <- rowSums(matrix_deconv_new_collapsed[,grepl("NPC", colnames(matrix_deconv_new_collapsed))])
matrix_deconv_new_collapsed <- matrix_deconv_new_collapsed[,!grepl("NPC", colnames(matrix_deconv_new_collapsed))]
matrix_deconv_new_collapsed$NPC <- npc_collapsed

matrix_deconv_new_collapsed_long <- as.data.frame(matrix_deconv_new_collapsed %>% pivot_longer(-c("Mixture"), names_to="CellType", values_to="Proportion"))


# go back to wide format
mat <- matrix_deconv_new %>%
  column_to_rownames("Mixture") %>%
  as.matrix()

mat_collapsed <- matrix_deconv_new_collapsed %>%
  column_to_rownames("Mixture") %>%
  as.matrix()

hc <- hclust(dist(mat), method = "ward.D2")

matrix_deconv_new_long$Mixture <- factor(
  matrix_deconv_new_long$Mixture,
  levels = rownames(mat)[hc$order]
)


getPalette = colorRampPalette(brewer.pal(7, "Set3"))
colourCount = length(unique(matrix_deconv_new_long$CellType))
colVec <- setNames(getPalette(colourCount),
                   sort(unique(matrix_deconv_new_long$CellType)))

matrix_deconv_new_long$CellType <- factor(matrix_deconv_new_long$CellType, levels=rev(names(colVec)))
matrix_deconv_new_long$Mixture <- factor(matrix_deconv_new_long$Mixture, levels=names(mat[,1])[order(mat[,1], decreasing=F)])


stackedBarplot <- ggplot(matrix_deconv_new_long, aes(fill=CellType, x=Mixture, y=Proportion))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 11, face="bold", hjust=0.5),
        legend.text  = element_text(size = 9),
        legend.key.size = unit(0.5, "lines"),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #legend.position="none",
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 10), ncol=1))+
  scale_fill_viridis_d(name="Cell types", option = "magma", direction=)+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("")+
  ylab("Percentage")+
  ggtitle("Predicted cell type proportions (CibersortX)")


pdf("deconvolution/plots/stackedBarplot_not_collapsed.pdf", width = 10, height = 8)
stackedBarplot
dev.off()

## collapse categories


matrix_deconv_new_collapsed_long$CellType <- factor(matrix_deconv_new_collapsed_long$CellType, levels=c("PNS","Pericytes","NPC","EN"))
matrix_deconv_new_collapsed_long$Mixture <- factor(matrix_deconv_new_collapsed_long$Mixture,
                                                   levels=names(mat_collapsed[,1])[order(mat_collapsed[,3], decreasing=T)])

matrix_deconv_new_collapsed_long$Condition <- metadata_bulk[match(matrix_deconv_new_collapsed_long$Mixture,
                                                                  metadata_bulk$seqsample),]$phenotype


vecCondition <- setNames(c("Control","ID"),
                         nm=c("con","in"))


matrix_deconv_new_collapsed_long$Condition2 <- unname(vecCondition[matrix_deconv_new_collapsed_long$Condition])


stackedBarplot2 <- ggplot(matrix_deconv_new_collapsed_long, aes(fill=CellType, x=Mixture, y=Proportion))+
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  facet_grid(~Condition2, scales = "free_x", space = "free_x")+
  theme(plot.title=element_text(size=18, face="bold", hjust=0.5),
        legend.title = element_text(size = 11, face="bold", hjust=0.5),
        legend.text  = element_text(size = 9),
        legend.key.size = unit(0.5, "lines"),
        axis.text.x = element_text(size=9, angle=90, vjust=0.5, hjust=1),
        #axis.text.x=element_blank(),
        #axis.ticks.x=element_blank(),
        #legend.position="none",
        axis.text.y = element_text(size=14),
        axis.title=element_text(size=12),
        strip.text.x = element_text(size = 10)) +
  guides(shape = guide_legend(override.aes = list(size = 10)),
         fill = guide_legend(override.aes = list(size = 10), ncol=1))+
  scale_fill_viridis_d(name="Cell types", option = "viridis", direction=)+
  scale_y_continuous(labels = scales::percent, breaks=seq(0,1,0.2))+
  xlab("")+
  ylab("Percentage")+
  ggtitle("Predicted cell type proportions (CibersortX)")

pdf("deconvolution/plots/figS3H.pdf", width = 12, height = 8)
stackedBarplot2
dev.off()


### edgeR/formal test

library(tidyverse)
library(scales)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)

#####
#####

pvalConverter <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.1)){
    symbolVec[vectorPval<0.1 & (!is.na(vectorPval))] <- "s"
  }
  
  if (any(vectorPval<0.05  & (!is.na(vectorPval)))){
    symbolVec[vectorPval<0.05 & (!is.na(vectorPval))] <- "*"
  }
  
  if (any(vectorPval<0.01 & (!is.na(vectorPval))) ){
    symbolVec[vectorPval<0.01 & (!is.na(vectorPval))] <- "**"
  }
  
  if (any(vectorPval<0.001 & (!is.na(vectorPval)))){
    symbolVec[vectorPval<0.001 & (!is.na(vectorPval))] <- "***"
  }
  
  if (any(is.na(vectorPval))){
    symbolVec[is.na(vectorPval)] <- "na"
  }
  
  
  return(symbolVec)
}

######
######

pvalConverter2 <- function(vectorPval){
  symbolVec <- rep("ns", length(vectorPval))
  
  if (any(vectorPval<0.1)){
    symbolVec[vectorPval<0.1 & (!is.na(vectorPval))] <- "p<0.1"
  }
  
  if (any(vectorPval<0.05  & (!is.na(vectorPval)))){
    symbolVec[vectorPval<0.05 & (!is.na(vectorPval))] <- "p<0.05"
  }
  
  if (any(vectorPval<0.01 & (!is.na(vectorPval))) ){
    symbolVec[vectorPval<0.01 & (!is.na(vectorPval))] <- "p<0.01"
  }
  
  if (any(vectorPval<0.001 & (!is.na(vectorPval)))){
    symbolVec[vectorPval<0.001 & (!is.na(vectorPval))] <- "p<0.001"
  }
  
  if (any(is.na(vectorPval))){
    symbolVec[is.na(vectorPval)] <- "na"
  }
  
  
  return(symbolVec)
}


### abbrev_not_collapsed

csSnapseed <- read.table("deconvolution/cibersortx/output/snapseed/abbreviation/CIBERSORTx_Adjusted.txt", header=T)
stopifnot(all(abs(rowSums(csSnapseed[,2:11])-1)<1e-8))


dataframeConverterCibersortX <- function(csSnapseed,  namesTo="CellTypes", valuesTo="Proportion", group="SetYourGroupColumn", group2="SetYourGroupColumn2"){
  
  to_rm <- (dim(csSnapseed)[2]-2):(dim(csSnapseed)[2])
  csSnapseed <- csSnapseed[,-to_rm]
  dfObj <- as.data.frame(csSnapseed %>% pivot_longer(-c(1), names_to=namesTo, values_to=valuesTo))
  colnames(dfObj)[1] <- "Samples"
  dfObj$group <- group
  dfObj$group2 <- group2
  return(dfObj)
  
}

csSnapseed_df <- dataframeConverterCibersortX(csSnapseed,  namesTo="CellTypes", valuesTo="Proportion", group="CibersortX", group2="Snapseed")
csSnapseed_df$CellType <- factor(csSnapseed_df$CellType, levels=rev(names(colVec)))

resultsDeconv <- ggplot(csSnapseed_df, aes(x=CellTypes, y=round(Proportion*100,2)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width = 0.2)+
  theme_bw()+
  ylab("Estimated percentage [%]")+
  xlab("")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14),
        plot.title=element_text(size=17, face="bold", color="black", hjust=0.5),
        legend.position="top",
        legend.title=element_text(size=15, color="black", face="bold"),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size=13))+
  ggtitle("Deconvolution (CibersortX)")



ggsave(filename = "deconvolution/plots/proportionsPerCellType_donor.pdf",
       plot = resultsDeconv,
       width = 7, height = 7)


######
######


library(SummarizedExperiment)
library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)
library(sva)
library(geneplotter)


deconvTab_matrix <- csSnapseed
rownames(deconvTab_matrix) <- deconvTab_matrix$Mixture
deconvTab_matrix$Mixture <- NULL
cols_to_rm <- (dim(deconvTab_matrix)[2]-2):(dim(deconvTab_matrix)[2])
deconvTab_matrix <- deconvTab_matrix[,-cols_to_rm]

metadata_bulk <- metadata_bulk[match(rownames(deconvTab_matrix),metadata_bulk$seqsample),]

## pass metadata
stopifnot(all(rownames(deconvTab_matrix)==metadata_bulk$seqsample))

# do I need to select only cell types present in a certain fraction?
#stopifnot(all(colMeans(deconvTab_matrix)>0.01))

table(metadata_bulk$sex, metadata_bulk$phenotype)
# con in
# f  10  5
# m   6  6



##Strategy 1 : edgeR/limma

library(limma)
library(sva)


tmp_deconvTab_matrix <- round(deconvTab_matrix[colMeans(deconvTab_matrix)>0.01]*10000)
  
metadata_bulk$phenotype <- factor(metadata_bulk$phenotype, levels=c("con","in"))
print(table(metadata_bulk$phenotype))

metadata_bulk$sex <- factor(metadata_bulk$sex, levels=c("f","m"))
print(table(metadata_bulk$phenotype))

metadata_bulk$rnu4.1 <- factor(metadata_bulk$rnu4.1, levels=c("N","Y"))
print(table(metadata_bulk$rnu4.1))

  
y.ab <- DGEList(t(as.matrix(tmp_deconvTab_matrix)), samples=metadata_bulk)
keep <- filterByExpr(y.ab, group=y.ab$samples$phenotype)
not_to_keep <- names(which(!keep))
  
if (length(not_to_keep)>0){
  tmp_deconvTab_matrix <- tmp_deconvTab_matrix[,-match(not_to_keep, colnames(tmp_deconvTab_matrix))]
}
  
mod <- model.matrix(~ phenotype + sex + astro_percentage + Neurog2 + rnu4.1, metadata_bulk)
mod0 <- model.matrix(~1 , metadata_bulk)

cor(model.matrix(~ phenotype + sex + astro_percentage + Neurog2 + rnu4.1, metadata_bulk)[, -1])
# phenotypein       sexm astro_percentage     Neurog2    rnu4.1Y
# phenotypein       1.00000000 0.16854997        0.2155665  0.06165274  0.6931818
# sexm              0.16854997 1.00000000        0.1602389  0.06683117  0.1685500
# astro_percentage  0.21556653 0.16023891        1.0000000 -0.37453983  0.3002171
# Neurog2           0.06165274 0.06683117       -0.3745398  1.00000000 -0.2894961
# rnu4.1Y           0.69318182 0.16854997        0.3002171 -0.28949615  1.0000000

#IQRs <- apply(assays(se.filt)$logCPM, 1, IQR)
# sv <- sva(t(as.matrix(tmp_deconvTab_matrix)),
#           mod=mod, mod0=mod0)
# if (sv$n.sv>0){
#   cnames <- c(colnames(mod), paste0("SV", 1:sv$n))
#   mod <- cbind(mod, sv$sv)
#   colnames(mod) <- cnames
#   head(mod, n=3)
# }
  
y.ab <- DGEList(t(as.matrix(tmp_deconvTab_matrix)), samples=metadata_bulk)
  
y.ab <- estimateDisp(y.ab, mod, robust=TRUE)
prdf <- cut(y.ab$prior.df, breaks=c(0, 1, 2, 3, 4, 5))
table(prdf, useNA="always")
  
v <- voom(y.ab, mod)
fit <- lmFit(v, mod)
fit <- eBayes(fit, robust=TRUE)
res <- decideTests(fit, p.value=0.1)

summary(res)
# (Intercept) phenotypein sexm astro_percentage Neurog2 rnu4.1Y
# Down             0           0    0                2       1       0
# NotSig           0           8    8                4       5       7
# Up               8           0    0                2       2       1

limmaRes_pheno <- topTable(fit, coef="phenotypein", n=Inf)
limmaRes_sex <- topTable(fit, coef="sexm", n=Inf)
limmaRes_astrop <- topTable(fit, coef="astro_percentage", n=Inf)
limmaRes_ngn2 <- topTable(fit, coef="Neurog2", n=Inf)
limmaRes_rnu41 <- topTable(fit, coef="rnu4.1Y", n=Inf)


vars_to_cbind <- ls()[grepl("limmaRes_", ls())]

allRes <- sapply(vars_to_cbind, function(x){
  
  tmp <- get(x)
  tmp_new <- tmp[,c("adj.P.Val","logFC")]
  tmp_new$CellType <- rownames(tmp_new)
  rownames(tmp_new) <- NULL
  tmp_new$ModelCoefficient <- gsub("limmaRes_","", x)
  
  return(tmp_new)
  
}, simplify=F)



allRes <- do.call("rbind", allRes); rownames(allRes) <- NULL

allRes$significant <- allRes$adj.P.Val<0.1


### volcano plot

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

colourCount = length(unique(allRes$CellType))
colVec <- setNames(getPalette(colourCount),
                   unique(allRes$CellType))

allRes$ModelCoefficient <- factor(
  allRes$ModelCoefficient,
  levels = c("astrop", "ngn2", "pheno", "sex","rnu41"),
  labels = c(
    "Astrocyte %",
    "Ngn2 Expression",
    "Phenotype (ID vs Ctrl)",
    "Sex (M vs F)",
    "RNU4.1 (YES vs NO)"
  )
)

plot_se_new_protocol <- ggplot(allRes, aes(x=logFC, y=-log10(adj.P.Val), col=CellType, shape=significant))+
  geom_point(size=4)+
  theme_bw()+
  ylab(expression(-log[10] * Pval)) +
  xlab("log2FC (per relative unit)")+
  ggtitle("Differential abundance (CibersortX): Values per coefficient")+
  facet_wrap(~ModelCoefficient, scales="free_x")+
  geom_text_repel(data=subset(allRes, significant==TRUE),
                  aes(x=logFC, y=-log10(adj.P.Val), label=CellType),
                  min.segment.length = Inf, seed = 123, box.padding = 0.3, col="black", size=3)+
  theme(axis.text.x=element_text(size=14, angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=17),
        plot.title=element_text(size=16, hjust=0.5),
        legend.text=element_text(size=14),  # Increase legend text size
        legend.title=element_text(size=15, face="bold"),  # Increase legend title size
        legend.key=element_rect(size=5),  # Increase the size of the legend key
        legend.key.size=unit(1.5, 'lines'),
        strip.text=element_text(size=15))+
  scale_color_manual(name="Cell types (Snapseed)", values=colVec)+
  guides(alpha = "none", text = "none")+
  geom_hline(yintercept=-log10(0.1) , linetype="dashed", col="grey")

ggsave(filename = "deconvolution/plots/volcanoPlot_coeff.pdf",
       plot = plot_se_new_protocol,
       width = 10, height = 10)


## heatmap

allRes$signifPlot <- pvalConverter2(allRes$adj.P.Val)


size_values=c("ns"=0,
              "p<0.1"=2.5,
              "p<0.05"=5.5,
              "p<0.01"=8)

myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
my_breaks <- seq(-4,10,4)


allRes$signifPlot <- factor(allRes$signifPlot, levels=c("ns","p<0.1","p<0.05","p<0.01"))


library(patchwork)

plots <- lapply(
  split(allRes, allRes$ModelCoefficient),
  function(df) {
    
    lim <- max(abs(df$logFC), na.rm = TRUE)
    
    ggplot(df, aes(x = CellType, y = 1, fill = logFC)) +
      geom_tile(color = "black") +
      geom_point(aes(size = signifPlot), color = "black") +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        limits = c(-lim, lim),
        name = "log2FC/rel.unit",
        midpoint = 0
      ) +
      scale_size_manual(values = size_values, name = "Adj.Pval") +
      labs(
        title = unique(df$ModelCoefficient),
        y = NULL,
        x = "Cell types"
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
)

global_plot <- wrap_plots(plots, ncol = 1)


ggsave(filename = "deconvolution/plots/log2FC_perCoeff.pdf",
       plot = global_plot,
       width = 12, height = 16)



### abbrev_collapsed

csSnapseed <- read.table("deconvolution/cibersortx/output/snapseed/abbreviation/CIBERSORTx_Adjusted.txt", header=T)
stopifnot(all(abs(rowSums(csSnapseed[,2:11])-1)<1e-8))


en_collapsed <- rowSums(csSnapseed[,grepl("EN", colnames(csSnapseed))])
matrix_deconv_new_collapsed <- csSnapseed[,!grepl("EN", colnames(csSnapseed))]
matrix_deconv_new_collapsed$EN <- en_collapsed

npc_collapsed <- rowSums(matrix_deconv_new_collapsed[,grepl("NPC", colnames(matrix_deconv_new_collapsed))])
matrix_deconv_new_collapsed <- matrix_deconv_new_collapsed[,!grepl("NPC", colnames(matrix_deconv_new_collapsed))]
matrix_deconv_new_collapsed$NPC <- npc_collapsed

matrix_deconv_new_collapsed <- matrix_deconv_new_collapsed[,c(1:3,7:8,4:6)]


dataframeConverterCibersortX <- function(csSnapseed,  namesTo="CellTypes", valuesTo="Proportion", group="SetYourGroupColumn", group2="SetYourGroupColumn2"){
  
  to_rm <- (dim(csSnapseed)[2]-2):(dim(csSnapseed)[2])
  csSnapseed <- csSnapseed[,-to_rm]
  dfObj <- as.data.frame(csSnapseed %>% pivot_longer(-c(1), names_to=namesTo, values_to=valuesTo))
  colnames(dfObj)[1] <- "Samples"
  dfObj$group <- group
  dfObj$group2 <- group2
  return(dfObj)
  
}

csSnapseed_df <- dataframeConverterCibersortX(matrix_deconv_new_collapsed,  namesTo="CellTypes", valuesTo="Proportion", group="CibersortX", group2="Snapseed")
csSnapseed_df$CellTypes <- factor(csSnapseed_df$CellTypes, levels=unique(csSnapseed_df$CellTypes))

resultsDeconv <- ggplot(csSnapseed_df, aes(x=CellTypes, y=round(Proportion*100,2)))+
  geom_boxplot(outlier.shape=NA)+
  geom_jitter(width = 0.2)+
  theme_bw()+
  ylab("Estimated percentage [%]")+
  xlab("")+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=7, color="black"),
        axis.text.y=element_text(size=12, color="black"),
        axis.title=element_text(size=14),
        plot.title=element_text(size=17, face="bold", color="black", hjust=0.5),
        legend.position="top",
        legend.title=element_text(size=15, color="black", face="bold"),
        legend.text=element_text(size=12),
        strip.text.x = element_text(size=13))+
  ggtitle("Deconvolution (CibersortX)")



ggsave(filename = "deconvolution/plots/proportionsPerCellType_donor_collapsed.pdf",
       plot = resultsDeconv,
       width = 7, height = 7)


######
######


library(SummarizedExperiment)
library(edgeR)
library(limma)
library(DESeq2)
library(ggplot2)
library(sva)
library(geneplotter)


#deconvTab_matrix <- csSnapseed
rownames(matrix_deconv_new_collapsed) <- matrix_deconv_new_collapsed$Mixture
matrix_deconv_new_collapsed$Mixture <- NULL
cols_to_rm <- (dim(matrix_deconv_new_collapsed)[2]-2):(dim(matrix_deconv_new_collapsed)[2])
matrix_deconv_new_collapsed <- matrix_deconv_new_collapsed[,-cols_to_rm]

metadata_bulk <- metadata_bulk[match(rownames(matrix_deconv_new_collapsed),metadata_bulk$seqsample),]

## pass metadata
stopifnot(all(rownames(matrix_deconv_new_collapsed)==metadata_bulk$seqsample))

# do I need to select only cell types present in a certain fraction?
#stopifnot(all(colMeans(matrix_deconv_new_collapsed)>0.01))

table(metadata_bulk$sex, metadata_bulk$phenotype)
# con in
# f  10  5
# m   6  6



##Strategy 1 : edgeR/limma

library(limma)
library(sva)

tmp_deconvTab_matrix <- round(matrix_deconv_new_collapsed[colMeans(matrix_deconv_new_collapsed)>0.01]*10000)

metadata_bulk$phenotype <- factor(metadata_bulk$phenotype, levels=c("con","in"))
print(table(metadata_bulk$phenotype))

metadata_bulk$sex <- factor(metadata_bulk$sex, levels=c("f","m"))
print(table(metadata_bulk$sex))

metadata_bulk$rnu4.1 <- factor(metadata_bulk$rnu4.1, levels=c("N","Y"))
print(table(metadata_bulk$rnu4.1))


y.ab <- DGEList(t(as.matrix(tmp_deconvTab_matrix)), samples=metadata_bulk)
keep <- filterByExpr(y.ab, group=y.ab$samples$phenotype)
not_to_keep <- names(which(!keep))

if (length(not_to_keep)>0){
  tmp_deconvTab_matrix <- tmp_deconvTab_matrix[,-match(not_to_keep, colnames(tmp_deconvTab_matrix))]
}

mod <- model.matrix(~ phenotype + sex + astro_percentage + Neurog2 + rnu4.1, metadata_bulk)
mod0 <- model.matrix(~1 , metadata_bulk)

cor(model.matrix(~ phenotype + sex + astro_percentage + Neurog2 + rnu4.1, metadata_bulk)[, -1])
#                   phenotypein     sexm  astro_percentage     Neurog2    rnu4.1Y
# phenotypein       1.00000000 0.16854997        0.2155665  0.06165274  0.6931818
# sexm              0.16854997 1.00000000        0.1602389  0.06683117  0.1685500
# astro_percentage  0.21556653 0.16023891        1.0000000 -0.37453983  0.3002171
# Neurog2           0.06165274 0.06683117       -0.3745398  1.00000000 -0.2894961
# rnu4.1Y           0.69318182 0.16854997        0.3002171 -0.28949615  1.0000000

#IQRs <- apply(assays(se.filt)$logCPM, 1, IQR)
# sv <- sva(t(as.matrix(tmp_deconvTab_matrix)),
#           mod=mod, mod0=mod0)
# if (sv$n.sv>0){
#   cnames <- c(colnames(mod), paste0("SV", 1:sv$n))
#   mod <- cbind(mod, sv$sv)
#   colnames(mod) <- cnames
#   head(mod, n=3)
# }

y.ab <- DGEList(t(as.matrix(tmp_deconvTab_matrix)), samples=metadata_bulk)

y.ab <- estimateDisp(y.ab, mod, robust=TRUE)
prdf <- cut(y.ab$prior.df, breaks=c(0, 1, 2, 3, 4, 5))
table(prdf, useNA="always")

v <- voom(y.ab, mod)
fit <- lmFit(v, mod)
fit <- eBayes(fit, robust=TRUE)
res <- decideTests(fit, p.value=0.1)

summary(res)
# (Intercept) phenotypein sexm astro_percentage Neurog2 rnu4.1Y
# Down             0           0    1                1       0       0
# NotSig           0           3    1                1       3       3
# Up               3           0    1                1       0       0

limmaRes_pheno <- topTable(fit, coef="phenotypein", n=Inf)
limmaRes_sex <- topTable(fit, coef="sexm", n=Inf)
limmaRes_astrop <- topTable(fit, coef="astro_percentage", n=Inf)
limmaRes_ngn2 <- topTable(fit, coef="Neurog2", n=Inf)
limmaRes_rnu41 <- topTable(fit, coef="rnu4.1Y", n=Inf)


vars_to_cbind <- ls()[grepl("limmaRes_", ls())]

allRes <- sapply(vars_to_cbind, function(x){
  
  tmp <- get(x)
  tmp_new <- tmp[,c("adj.P.Val","logFC")]
  tmp_new$CellType <- rownames(tmp_new)
  rownames(tmp_new) <- NULL
  tmp_new$ModelCoefficient <- gsub("limmaRes_","", x)
  
  return(tmp_new)
  
}, simplify=F)



allRes <- do.call("rbind", allRes); rownames(allRes) <- NULL

allRes$significant <- allRes$adj.P.Val<0.1


### volcano plot

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

colourCount = length(unique(allRes$CellType))
colVec <- setNames(getPalette(colourCount),
                   unique(allRes$CellType))

allRes$ModelCoefficient <- factor(
  allRes$ModelCoefficient,
  levels = c("astrop", "ngn2", "pheno", "sex","rnu41"),
  labels = c(
    "Astrocyte %",
    "Ngn2 Expression",
    "Phenotype (ID vs Ctrl)",
    "Sex (M vs F)",
    "RNU4.1 (YES vs NO)"
  )
)

plot_se_new_protocol <- ggplot(allRes, aes(x=logFC, y=-log10(adj.P.Val), col=CellType, shape=significant))+
  geom_point(size=4)+
  theme_bw()+
  ylab(expression(-log[10] * Pval)) +
  xlab("log2FC (per relative unit)")+
  ggtitle("Differential abundance (CibersortX): Values per coefficient")+
  facet_wrap(~ModelCoefficient, scales="free_x")+
  geom_text_repel(data=subset(allRes, significant==TRUE),
                  aes(x=logFC, y=-log10(adj.P.Val), label=CellType),
                  min.segment.length = Inf, seed = 123, box.padding = 0.3, col="black", size=6)+
  theme(axis.text.x=element_text(size=14, angle=90, hjust=0.5, vjust=0.5),
        axis.text.y=element_text(size=14),
        axis.title=element_text(size=17),
        plot.title=element_text(size=16, hjust=0.5),
        legend.text=element_text(size=14),  # Increase legend text size
        legend.title=element_text(size=15, face="bold"),  # Increase legend title size
        legend.key=element_rect(size=5),  # Increase the size of the legend key
        legend.key.size=unit(1.5, 'lines'),
        strip.text=element_text(size=15))+
  scale_color_manual(name="Cell types (Snapseed)", values=colVec)+
  guides(alpha = "none", text = "none")+
  geom_hline(yintercept=-log10(0.1) , linetype="dashed", col="grey")

ggsave(filename = "deconvolution/plots/volcanoPlot_coeff_collapsed.pdf",
       plot = plot_se_new_protocol,
       width = 10, height = 10)


## heatmap

allRes$signifPlot <- pvalConverter2(allRes$adj.P.Val)


size_values=c("ns"=0,
              "p<0.1"=2,
              "p<0.05"=4,
              "p<0.01"=6,
              "p<0.001"=8)



myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
my_breaks <- seq(-4,10,4)


allRes$signifPlot <- factor(allRes$signifPlot, levels=c("ns","p<0.1","p<0.05","p<0.01","p<0.001"))

library(patchwork)

plots <- lapply(
  split(allRes, allRes$ModelCoefficient),
  function(df) {
    
    lim <- max(abs(df$logFC), na.rm = TRUE)
    
    ggplot(df, aes(x = CellType, y = 1, fill = logFC)) +
      geom_tile(color = "black") +
      geom_point(aes(size = signifPlot), color = "black") +
      scale_fill_gradient2(
        low = "blue", mid = "white", high = "red",
        limits = c(-lim, lim),
        name = "log2FC/rel.unit",
        midpoint = 0
      ) +
      scale_size_manual(values = size_values, name = "Adj.Pval") +
      labs(
        title = unique(df$ModelCoefficient),
        y = NULL,
        x = "Cell types"
      ) +
      theme_bw() +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
)

global_plot <- wrap_plots(plots, ncol = 1)


ggsave(filename = "deconvolution/plots/log2FC_perCoeff_collapsed.pdf",
       plot = global_plot,
       width =4, height = 14)



## check at significant hits

df_val <- matrix_deconv_new_collapsed[,-2]
df_val$sample <- rownames(df_val)
rownames(df_val) <- NULL

df_val$sex <- metadata_bulk[match(df_val$sample, metadata_bulk$seqsample),]$sex
df_val$astro_percentage <- metadata_bulk[match(df_val$sample, metadata_bulk$seqsample),]$astro_percentage
df_val$phenotype <- metadata_bulk[match(df_val$sample, metadata_bulk$seqsample),]$phenotype

df_val_long <- as.data.frame(df_val %>% pivot_longer(-c("sample","sex","astro_percentage","phenotype"), names_to="CellTypes", values_to="Proportion"))

library(ggpubr)

sex_plot <- ggplot(df_val_long, aes(x=sex, y=Proportion))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~CellTypes)+
  theme_bw()+
  ggtitle("Sex differences (per cell type)")+
  xlab("Sex")+
  ylab("Inferred cell type proportion")+
  stat_compare_means()

ggsave(filename = "deconvolution/plots/validation_sex_plot_collapsed.pdf",
       plot = sex_plot,
       width =6, height = 4.5)

condition_plot <- ggplot(df_val_long, aes(x=phenotype, y=Proportion))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~CellTypes)+
  theme_bw()+
  ggtitle("Condition differences (per cell type)")+
  xlab("Condition")+
  ylab("Inferred cell type proportion")+
  stat_compare_means()

ggsave(filename = "deconvolution/plots/figS3I.pdf",
       plot = condition_plot,
       width =6, height = 4.5)



astro_perc <- ggplot(df_val_long, aes(x = astro_percentage, y = Proportion)) +
  geom_point(pch=21, size=3) +
  facet_wrap(~CellTypes) +
  theme_bw() +
  geom_smooth(method = "lm", se = FALSE) +
  stat_cor(method = "pearson",
           label.x.npc = "left",
           label.y.npc = 1,
           size = 3) +
  xlab("Astrocyte cells [%]") +
  ylab("Inferred cell type proportion") +
  ggtitle("Linear effects of astrocyte abundance on cell type composition")


ggsave(filename = "deconvolution/plots/figS3J.pdf",
       plot = astro_perc,
       width =8, height = 5.5)









